import numpy as np
from scipy.sparse.linalg import cg
from scipy.sparse import coo_matrix

import cedar


class HeatTransfer(cedar.Model):
    """
    Heat conduction model using FVM.

    Solves the heat transfer equation on a 3D mesh with region-dependent
    materials, heat sources, and boundary conditions.

    Parameters
    ----------
    mesh : cedar.Mesh3D
        The mesh that the heat conduction equation will be solved on.

    Attributes
    ----------
    mesh : cedar.Mesh3D
        The mesh that the heat conduction equation will be solved on.
    materials : dict[str, cedar.Material]
        Material properties objects. One entry must be provided for each region.
        You can add entries with HeatTransfer.add_material().
    bcs : dict[str, cedar.HeatTransferBC]
        Boundary Conditions. The keys are boundary names in HeatTransfer.mesh.
        Any boundary without an entry listed here will default to
        cedar.AdiabaticBC().
    sources: dict[str, cedar.HeatTransferSource]
        Source Terms. The keys are region names in HeatTransfer.mesh. Any region
        without an entry will have zero heat generation.
    T : cedar.Field
        Temperature in [K].
    Qgen : cedar.Field
        Heat generation in [W].
    volQgen : cedar.Field
        Volumetric heat generation in [W/m^3].
    k : cedar.Field
        Thermal conductivity in [W/m-K]
    rho : cedar.Field
        Density in [kg/m3]
    cp : cedar.Field
        Specific heat capacity in [J/kg-K].
    J : cedar.Field
        Boundary heat flux in [W/m^2].
    Qdot : cedar.Field
        Boundary heat transfer in [W].

    """

    def __init__(
        self,
        mesh: cedar.Mesh3D = None,
    ):
        self.mesh: cedar.Mesh3D = mesh
        self.materials: dict[str, cedar.Material] = {}
        self.bcs: dict[str, cedar.HeatTransferBC] = {}
        self.sources: dict[str, cedar.HeatTransferSource] = {}

        self.T: cedar.Field = self._add_field(
            "T", "HeatTransfer", 300, is_on_boundaries=True
        )

        self.Qgen: cedar.Field = self._add_field(
            "Qgen", "HeatTransfer", 0, is_on_boundaries=True
        )
        self.volQgen: cedar.Field = self._add_field("volQgen", "HeatTransfer", 0)
        self.k: cedar.Field = self._add_field("k", "HeatTransfer", 1)
        self.rho: cedar.Field = self._add_field("rho", "HeatTransfer", 1)
        self.cp: cedar.Field = self._add_field("cp", "HeatTransfer", 1)

        self.J: cedar.Field = self._add_field(
            "J", "HeatTransfer", 0, is_on_regions=False, is_on_boundaries=True
        )
        self.Qdot: cedar.Field = self._add_field(
            "Qdot", "HeatTransfer", 0, is_on_regions=False, is_on_boundaries=True
        )

    def add_source(self, region: str, source: cedar.HeatTransferSource):
        """
        Add a heat generation term to the heat transfer model.

        Parameters
        ----------
        region : str
            Region to apply heat generation term.
        source : cedar.HeatTransferSource
            Heat generation source.
        """

        self.sources[region] = source

    def set_bc(self, boundary: str, bc: cedar.HeatTransferBC):
        """
        Set boundary condition of the heat transfer model.

        Parameters
        ----------
        boundary : str
            Boundary to apply boundary condition.
        bc : cedar.HeatTransferBC
            Boundary condition.

        """
        self.bcs[boundary] = bc

    def set_material(self, region: str, material: cedar.Material):
        """
        Set material in the heat transfer model.

        Parameters
        ----------
        region : str
            Region to set material.
        material : cedar.Material
            Material properties object.
        """
        self.materials[region] = material

    def _A_b(self, k_face, Qgen, cp, T_ic, dt):
        N = self.mesh.N_cells
        b = Qgen

        rows = []
        cols = []
        data = []

        if dt is not None:
            C = self._mass * cp / dt

            rows.extend(range(N))
            cols.extend(range(N))
            data.extend(C)

            b += T_ic * C

        # -----------------------------------------
        #  INTERIOR FACES
        # -----------------------------------------
        face_data = zip(
            self.mesh.face_cells_i[self.mesh.face_is_interior],
            k_face[self.mesh.face_is_interior],
            self._area_over_d[self.mesh.face_is_interior],
        )

        for cell_ids, k, area_over_d in face_data:
            C = k * area_over_d
            c1, c2 = cell_ids

            # A[c1][c1] += C
            rows.append(c1)
            cols.append(c1)
            data.append(C)
            # A[c1][c2] += -C
            rows.append(c1)
            cols.append(c2)
            data.append(-C)
            # A[c2][c2] += C
            rows.append(c2)
            cols.append(c2)
            data.append(C)
            # A[c2][c1] += -C
            rows.append(c2)
            cols.append(c1)
            data.append(-C)

        # -----------------------------------------
        #  BOUNDARY FACES
        # -----------------------------------------
        for boundary in self.mesh.boundaries:
            boundary_i = self.mesh.boundary_i[boundary]

            T = self.T.get(boundary)
            k = k_face[boundary_i]
            d = self.mesh.face_d[boundary_i, 0]
            area = self.mesh.face_areas[boundary_i]

            c1 = self.mesh.face_cells_i[boundary_i, 0]

            LHS, RHS = self.bcs[boundary]._ht_contribution(T, k, d, area)

            rows.extend(c1)
            cols.extend(c1)
            data.extend(LHS)
            b[c1] += RHS

        A = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

        return A, b

    def _cp_by_cell(self, T: np.ndarray):
        """
        Compute cell-wise specific heat capacity.
        """
        cp = np.zeros(self.mesh.N_cells, dtype=np.float64)

        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            cp[mask] = mat.cp(T[mask])

        return cp

    def _build(self):
        for boundary in self.mesh.boundaries:
            if boundary not in self.bcs:
                self.bcs[boundary] = cedar.AdiabaticBC()

        # We don't need to add a "default" heat source because it defaults to
        # zero unless populated by an existing source.

    def _initialize(self):
        self._area_over_d: np.ndarray = self.mesh.face_areas / np.sum(
            self.mesh.face_d, axis=1
        )

        self._rho = self._rho_by_cell()
        self._mass = self._rho * self.mesh.cell_vols

    def _iterate(self, res_reduc: float, dt: float) -> float:
        T = self.T.get()
        T_ic = self.T._get_ic()

        Qgen = self._Qgen_from_sources()

        volQgen = Qgen / self.mesh.cell_vols

        k_cell = self._k_by_cell(T)
        k_face = self._k_by_face(k_cell)
        cp = self._cp_by_cell(T)

        A, b = self._A_b(k_face, Qgen, cp, T_ic, dt)

        res = cedar.functions.residual(A, b, T)

        tgt_res = res_reduc * res

        T, _ = cg(A, b, T, rtol=tgt_res)

        self.T._set(T)
        self.Qgen._set(Qgen)
        self.volQgen._set(volQgen)
        self.k._set(k_cell)
        self.rho._set(self._rho)
        self.cp._set(cp)

        for boundary in self.mesh.boundaries:
            face_i = self.mesh.boundary_i[boundary]

            T_cell = T[self.mesh.face_cells_i[face_i, 0]]
            d = self.mesh.face_d[face_i, 0]

            T_boundary_face = self.bcs[boundary]._boundary_T(T_cell, k_face[face_i], d)
            self.T._set(T_boundary_face, boundary)

            J_boundary_face = self.bcs[boundary]._boundary_J(T_cell, k_face[face_i], d)
            self.J._set(J_boundary_face, boundary)

            Qdot_boundary_face = J_boundary_face * self.mesh.face_areas[face_i]
            self.Qdot._set(Qdot_boundary_face, boundary)

        return res

    def _k_by_cell(self, T: np.ndarray):
        """
        Compute cell-wise thermal conductivity.
        """
        k = np.empty(self.mesh.N_cells)

        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            k[mask] = mat.k(T[mask])

        return k

    def _k_by_face(self, k_cell):
        """
        Compute face thermal conductivity from adjacent cells.
        """
        k_face = np.empty(self.mesh.N_faces, dtype=np.float64)

        f = self.mesh.face_is_interior
        c1 = self.mesh.face_cells_i[f, 0]
        c2 = self.mesh.face_cells_i[f, 1]
        w = self.mesh.face_w[f]

        k_face[f] = self._k_face_interior(k_cell[c1], k_cell[c2], w)

        # boundary faces
        bf = ~f
        k_face[bf] = k_cell[self.mesh.face_cells_i[bf, 0]]

        return k_face

    def _k_face_interior(
        self, k1: np.ndarray, k2: np.ndarray, w: np.ndarray
    ) -> np.ndarray:
        return k1 * k2 / (k2 * (1 - w) + k1 * w)

    def _Qgen_from_sources(self) -> np.ndarray:
        Qgen = np.zeros(self.mesh.N_cells, dtype=np.float64)

        for region, source in self.sources.items():
            region_i = self.mesh.region_i[region]
            cell_vols = self.mesh.cell_vols[region_i]

            Qgen[region_i] = source._ht_contribution(cell_vols)

        return Qgen

    def _rho_by_cell(self):
        """
        Compute cell-wise density.

        :return: Density per cell.
        """
        rho = np.zeros(self.mesh.N_cells, dtype=np.float64)

        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            rho[mask] = mat.rho_rt()

        return rho
