from abc import ABC, abstractmethod
from typing import Callable
import numpy as np
from scipy.sparse.linalg import cg
from scipy.sparse import coo_matrix

import cedar


class HeatTransferBC(cedar.BC, ABC):
    """
    Abstract base class for :class:`cedar.HeatTransfer` boundary conditions.

    Defines the interface required for coupling boundary conditions into
    the heat transfer solver.
    """
    def __init__(self):
        """
        Initialize the boundary condition base class.
        """
        self.boundary: str = None
        self.mesh: cedar.Mesh = None

    @abstractmethod
    def boundary_B(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Compute the heat transfer through boundary.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Heat transfer through boundary.
        """
        pass

    @abstractmethod
    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Compute the boundary temperature values.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Boundary temperature values.
        """
        pass

    def initialize(self, t_start: float):
        """
        Initialize mesh-dependent boundary quantities.

        :param t_start: Simulation start time.
        """
        face_idx = self.mesh.boundary_i[self.boundary]

        self.area_over_d = (
            self.mesh.face_areas[face_idx] /
            self.mesh.face_d[face_idx, 0]
        ).astype(np.float64)

        self.update(t_start)

    @abstractmethod
    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Return matrix contributions for A, b assembly.

        :param T: Boundary temperature values.
        :param k: Face thermal conductivity.
        :return: (LHS, RHS) diagonal and source term contributions.
        """
        pass

    @abstractmethod
    def update(self, t: float):
        """
        Update the boundary condition.

        :param t: Current simulation time.
        """
        pass

    def evaluate_callable(self, bc_callable: Callable, t: float):
        """
        Evaluate spatially and temporally varying temperature function.

        :param t: Current simulation time.
        :return: Boundary temperature values.
        """
        values = np.zeros(self.mesh.boundary_N[self.boundary], dtype = np.float64)

        for i, face_i in enumerate(self.mesh.boundary_i[self.boundary]):
            x = self.mesh.face_centers[face_i][0]
            y = self.mesh.face_centers[face_i][1]
            z = self.mesh.face_centers[face_i][2]

            values[i] = bc_callable(x, y, z, t)

        return values

class AdiabaticBC(HeatTransferBC):
    """
    Adiabatic (zero-flux) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces zero normal heat flux across the boundary.
    """
    def __init__(self):
        """
        Initialize the adiabatic boundary condition.
        """
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

    def boundary_B(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Compute the heat transfer through boundary.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Heat transfer through boundary.
        """
        return np.zeros_like(T_cell)

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Temperature at boundary.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Zero source term contribution.
        """
        return np.copy(T_cell)

    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Return matrix contributions for A, b assembly.

        :param T: Boundary temperature values in [K].
        :param k: Face thermal conductivity.
        :return: (LHS, RHS) diagonal and source term contributions.
        """

        LHS = np.zeros(self.mesh.boundary_N[self.boundary], dtype=np.float64)
        RHS = np.zeros(self.mesh.boundary_N[self.boundary], dtype=np.float64)

        return LHS, RHS

    def update(self, t: float):
        """
        Update the boundary temperature for the current time.

        :param t: Current simulation time.
        """
        pass

class FixedTemperatureBC(HeatTransferBC):
    """
    Fixed-temperature (Dirichlet) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces a prescribed boundary temperature, optionally varying in space
    and time.
    """
    def __init__(self, T: float | Callable):
        """
        Initialize the fixed temperature boundary condition.

        :param T: Boundary temperature value or callable T(x, y, z, t).
        """
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

        self.T: float | Callable = T

        self.value: np.ndarray = None
        self.area_over_d: np.ndarray = None

    def boundary_B(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Compute the heat transfer through boundary.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Heat transfer through boundary.
        """
        return k*self.area_over_d*(T_cell-self.value)

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Temperature at boundary.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Zero source term contribution.
        """
        return self.value

    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Return matrix contributions for A, b assembly.

        :param T: Boundary temperature values in [K].
        :param k: Face thermal conductivity.
        :return: (LHS, RHS) diagonal and source term contributions.
        """
        LHS = k*self.area_over_d
        RHS = k*self.area_over_d * self.value

        return LHS, RHS

    def update(self, t: float):
        """
        Update the boundary temperature for the current time.

        :param t: Current simulation time.
        """
        if callable(self.T):
            self.value = self.evaluate_callable(self.T, t)

        else:
            self.value = np.full(self.mesh.boundary_N[self.boundary], self.T)
    
class FixedFluxBC(HeatTransferBC):
    """
    Fixed-flux (Neumann) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces a prescribed boundary flux, optionally varying in space
    and time.
    """
    def __init__(self, J: float | Callable):
        """
        Initialize the fixed-flux boundary condition.

        :param J: Boundary flux value or callable J(x, y, z, t) in [W/m^2].
        """
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

        self.J: float | Callable = J

        self.value: np.ndarray = None
        self.area_over_d: np.ndarray = None

    def boundary_B(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Compute the heat transfer through boundary.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Heat transfer through boundary.
        """
        face_i = self.mesh.boundary_i[self.boundary]
        return self.value * self.mesh.face_areas[face_i]

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Temperature at boundary.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Zero source term contribution.
        """
        face_i = self.mesh.boundary_i[self.boundary]
        return T_cell - self.mesh.face_d[face_i, 0]*self.value/k

    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Return matrix contributions for A, b assembly.

        :param T: Boundary temperature values in [K].
        :param k: Face thermal conductivity.
        :return: (LHS, RHS) diagonal and source term contributions.
        """
        face_i = self.mesh.boundary_i[self.boundary]
        LHS = np.zeros_like(T)
        RHS = -self.value*self.mesh.face_areas[face_i]

        return LHS, RHS

    def update(self, t: float):
        """
        Update the boundary temperature for the current time.

        :param t: Current simulation time.
        """
        if callable(self.J):
            self.value = self.evaluate_callable(self.J, t)

        else:
            self.value = np.full(self.mesh.boundary_N[self.boundary], self.J)

class HeatSource(cedar.Source):
    """
    Source for :class:`cedar.HeatTransfer`

    Distributes a total or volumetric heat input over the mesh using an
    optional spatial shape function and optional time dependence.
    """
    def __init__(
            self,
            Qdot: float | dict[str, float] | Callable = None,
            volQdot: float | dict[str, float] | Callable = None,
            shape_func: Callable = None
        ):
        """
        Initialize the heat source.

        :param Qdot: Total heat input (global, per-region, or time-dependent).
        :param volQdot: Volumetric heat input (global, per-region, or time-dependent).
        :param shape_func: Spatial weighting function f(x, y, z) or per-region mapping.
        """
        self.Qdot = Qdot
        self.volQdot = volQdot
        self.shape_func = shape_func

        self.mesh: cedar.Mesh3D = None

        self.Qdot_func: Callable = None
        self._value: np.ndarray = None
        self.w: np.ndarray = None

    def initialize(self, t_start: float):
        """
        Initialize the source distribution and evaluate initial value.

        :param t_start: Simulation start time.
        """
        # Spatial weighting (normalized)
        self.w = self._get_weighting()

        # Allocate cell-wise heat source
        self._value = np.zeros(self.mesh.N_cells, dtype = np.float64)

        # Convert volumetric source to total source if needed
        if self.volQdot is not None:
            if isinstance(self.volQdot, dict):
                self.Qdot = {}

                for region in self.mesh.regions:
                    self.Qdot[region] = self.volQdot[region] * self.mesh.region_vol[region]

            elif callable(self.volQdot):
                self.Qdot_func = lambda t: sum(self.mesh.region_vol.values())*self.volQdot(t)

            else:
                self.Qdot = self.volQdot * sum(self.mesh.region_vol.values())

        # Region-wise total heat input
        if isinstance(self.Qdot, dict):
            for region in self.mesh.regions:
                if region not in self.Qdot:
                    self.Qdot[region] = 0

                region_i = self.mesh.region_i[region]

                w = self.w[region_i]
                w = w / np.sum(w)

                self._value[region_i] = self.Qdot[region] * w

        # Time-dependent total heat input
        elif callable(self.Qdot):
            self.Qdot_func = self.Qdot

        # Static global total heat input
        else:
            self._value = self.Qdot * self.w

        self.step(t_start)

    def get_Qdot(self):
        """
        Return the current cell-wise heat source.

        :return: Heat source per cell.
        """
        return np.copy(self._value)

    def step(self, t: float):
        """
        Update the source value for the current time.

        :param t: Current simulation time.
        """
        if self.Qdot_func is not None:
            self._value = self.Qdot_func(t)*self.w

    def _get_weighting(self):
        """
        Compute normalized spatial weighting over the mesh.

        :return: Cell-wise weighting array summing to one.
        """
        # CALCULATE SHAPE
        if self.shape_func is None:
            w = np.ones(self.mesh.N_cells, dtype = np.float64)/self.mesh.N_cells

        else:
            w = np.zeros(self.mesh.N_cells, dtype = np.float64)

            for region in self.mesh.regions:
                region_i = self.mesh.region_i[region]

                if isinstance(self.shape_func, dict):
                    if region in self.shape_func:
                        shape_func = self.shape_func[region]

                    else:
                        shape_func = lambda x, y, z: 1

                else:
                    shape_func = self.shape_func

                for i in region_i:
                    x = self.mesh.cell_centers[i][0]
                    y = self.mesh.cell_centers[i][1]
                    z = self.mesh.cell_centers[i][2]

                    w[i] = shape_func(x, y, z)

        w = w * self.mesh.cell_vols

        return w / np.sum(w)
        
class HeatTransfer(cedar.Model):
    """
    Heat conduction model using FVM.

    Solves the heat equation on a 3D mesh with region-dependent materials, heat
    sources, and boundary conditions.
    """

    def __init__(
            self,
            mesh: cedar.Mesh3D = None,
            materials: dict[str, cedar.Material] = {},
            bc: dict[str, cedar.BC] = {},
            source: HeatSource = HeatSource(Qdot = 0)
    ):
        """
        Initialize the heat transfer model.

        :param mesh: Computational mesh.
        :param materials: Mapping of region name to material model.
        :param bc: Mapping of boundary name to boundary condition.
        :param source: Volumetric heat source model.
        """
        self.mesh: cedar.Mesh3D = mesh
        self.materials: dict[str, cedar.Material] = materials
        self.bc: dict[str, cedar.HeatTransferBC] = bc
        self.source: HeatSource = source

        self.T: cedar.Field = self.add_field("T", requires_ic = True, is_on_boundaries = True)
        self.Qdot: cedar.Field = self.add_field("Qdot")
        self.volQdot: cedar.Field = self.add_field("volQdot")
        self.k: cedar.Field = self.add_field("k")
        self.rho: cedar.Field = self.add_field("rho")
        self.cp: cedar.Field = self.add_field("cp")

        self.B: cedar.Field = self.add_field("B", is_on_regions = False, is_on_boundaries = True)

    def cp_by_cell(self, T: np.ndarray):
        """
        Compute cell-wise specific heat capacity.

        :param T: Cell temperature array.
        :return: Specific heat per cell.
        """
        cp = np.empty(self.mesh.N_cells)

        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            cp[mask] = mat.cp(T[mask])

        return cp

    def initialize(self, t_start: float):
        """
        Perform one-time initialization before solve.

        :param t_start: Simulation start time.
        """
        self._area_over_L: np.ndarray = self.mesh.face_areas / np.sum(self.mesh.face_d, axis=1)

        for boundary in self.mesh.boundaries:
            if boundary not in self.bc:
                self.bc[boundary] = AdiabaticBC()

        self._rho = self.rho_by_cell()
        self._mass = self._rho*self.mesh.cell_vols

    def iterate(self, res_reduc: float, dt: float) -> float:
        """
        Advance the solution by one nonlinear iteration.

        :param res_reduc: Target relative residual reduction.
        :param dt: Time step size.
        :return: Residual norm.
        """
        T = self.T.as_continuous_cell_value()

        Qdot = self.source.get_Qdot()
        volQdot = Qdot / self.mesh.cell_vols

        k_cell = self.k_by_cell(T)
        k_face = self.k_by_face(k_cell)
        cp = self.cp_by_cell(T)

        A, b = self._A_b(k_face, Qdot, cp, dt)

        res = cedar.functions.residual(A, b, T)
        
        tgt_res = res_reduc * res

        T, _ = cg(A, b, T, rtol = tgt_res)

        for region in self.mesh.regions:
            region_i = self.mesh.region_i[region]

            self.T.values[region][:] = T[region_i]
            self.Qdot.values[region][:] = Qdot[region_i]
            self.volQdot.values[region][:] = volQdot[region_i]
            self.k.values[region][:] = k_cell[region_i]
            self.rho.values[region][:] = self._rho[region_i]
            self.cp.values[region][:] = cp[region_i]

        for boundary in self.mesh.boundaries:
            face_i = self.mesh.boundary_i[boundary]
            T_cell = T[self.mesh.face_cells_i[face_i, 0]]

            self.T.values[boundary][:] = self.bc[boundary].boundary_T(T_cell, k_face[face_i])
            self.B.values[boundary][:] = self.bc[boundary].boundary_B(T_cell, k_face[face_i])

        return res

    def k_by_cell(self, T: np.ndarray):
        """
        Compute cell-wise thermal conductivity.

        :param T: Cell temperature array.
        :return: Thermal conductivity per cell.
        """
        k = np.empty(self.mesh.N_cells)
        
        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            k[mask] = mat.k(T[mask])

        return k

    def k_by_face(self, k_cell):
        """
        Compute face thermal conductivity from adjacent cells.

        Interior faces use harmonic averaging; boundary faces use the
        adjacent cell value.

        :param k_cell: Cell thermal conductivity.
        :return: Face thermal conductivity.
        """
        k_face = np.empty(self.mesh.N_faces, dtype = np.float64)

        f = self.mesh.face_is_interior
        c1 = self.mesh.face_cells_i[f, 0]
        c2 = self.mesh.face_cells_i[f, 1]
        w = self.mesh.face_w[f]

        k_face[f] = self._k_face_interior(k_cell[c1], k_cell[c2], w)

        # boundary faces
        bf = ~f
        k_face[bf] = k_cell[self.mesh.face_cells_i[bf, 0]]

        return k_face

    def rho_by_cell(self):
        """
        Compute cell-wise density.

        :return: Density per cell.
        """
        rho = np.empty(self.mesh.N_cells)
        
        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            rho[mask] = mat.rho_rt()

        return rho

    def _A_b(self, k_face, Qdot, cp, dt):
        """
        Assemble the linear system A T = b.

        :param k_face: Face thermal conductivity.
        :param Qdot: Cell heat source.
        :param cp: Cell specific heat capacity.
        :param dt: Time step size.
        :return: Sparse matrix A and RHS vector b.
        """
        N = self.mesh.N_cells
        b = Qdot

        rows = []
        cols = []
        data = []

        if dt is not None:
            C = self._mass * cp / dt
            T_ic = self.T.ic_as_continuous_cell_value()

            rows.extend(range(N)); cols.extend(range(N)); data.extend(C)

            b += T_ic * C

        # -----------------------------------------
        #  INTERIOR FACES
        # -----------------------------------------
        face_data = zip(
            self.mesh.face_cells_i[self.mesh.face_is_interior],
            k_face[self.mesh.face_is_interior],
            self._area_over_L[self.mesh.face_is_interior]
        )

        for cell_ids, k, area_over_L in face_data:
            C = k * area_over_L
            c1, c2 = cell_ids

            # A[c1][c1] += C
            rows.append(c1); cols.append(c1); data.append( C)
            # A[c1][c2] += -C
            rows.append(c1); cols.append(c2); data.append(-C)
            # A[c2][c2] += C
            rows.append(c2); cols.append(c2); data.append( C)
            # A[c2][c1] += -C
            rows.append(c2); cols.append(c1); data.append(-C)

        # -----------------------------------------
        #  BOUNDARY FACES
        # -----------------------------------------
        for boundary in self.mesh.boundaries:
            boundary_i = self.mesh.boundary_i[boundary]
            T_boundary = self.T.values[boundary]
            k_boundary = k_face[boundary_i]

            c1 = self.mesh.face_cells_i[boundary_i, 0]

            LHS, RHS = self.bc[boundary].matrix_contribution(T_boundary, k_boundary)

            rows.extend(c1); cols.extend(c1); data.extend(LHS)
            b[c1] += RHS

        A = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

        return A, b

    def _k_face_interior(self, k1, k2, w):
        return 1/((1-w)/k1 + w/k2)