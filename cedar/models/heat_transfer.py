from abc import ABC, abstractmethod
from typing import Callable
import numpy as np
from scipy.sparse.linalg import cg
from scipy.sparse import coo_matrix

import cedar


class HeatTransferBC(cedar.BC, ABC):
    def __init__(self):
        self.boundary: str = None
        self.mesh: cedar.Mesh = None

    @abstractmethod
    def boundary_B(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def initialize(self, t_start: float):
        pass

    @abstractmethod
    def iterate(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        pass

    @abstractmethod
    def step(self, t: float):
        pass

class AdiabaticBC(HeatTransferBC):
    def __init__(self):
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

    def boundary_B(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return np.zeros_like(T_cell)

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return T_cell

    def initialize(self, t_start: float):
        pass

    def iterate(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        LHS = np.zeros(self.mesh.boundary_N[self.boundary], dtype=np.float64)
        RHS = np.zeros(self.mesh.boundary_N[self.boundary], dtype=np.float64)

        return LHS, RHS

    def step(self, t: float):
        pass

class FixedTemperatureBC(HeatTransferBC):
    def __init__(self, T: float | Callable):
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

        self.T: float | Callable = T

        self._value: np.ndarray = None
        self._area_over_L: np.ndarray = None

    def boundary_B(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return k*self._area_over_L*(T_cell-self._value)

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return self._value

    def initialize(self, t_start: float):
        self._area_over_L = np.zeros(self.mesh.boundary_N[self.boundary], dtype=np.float64)

        for i, face_i in enumerate(self.mesh.boundary_i[self.boundary]):
            self._area_over_L[i] = self.mesh.face_areas[face_i] / self.mesh.face_d[face_i][0]

        self.step(t_start)

    def iterate(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        LHS = k*self._area_over_L
        RHS = k*self._area_over_L * self._value

        return LHS, RHS

    def step(self, t: float):
        if callable(self.T):
            self._value = self._evaluate_func(t)

        else:
            self._value = np.full(self.mesh.boundary_N[self.boundary], self.T)
    
    def _evaluate_func(self, t: float):
        values = np.zeros(self.mesh.boundary_N[self.boundary], dtype = np.float64)

        for i, face_i in enumerate(self.mesh.boundary_i[self.boundary]):
            x = self.mesh.face_centers[face_i][0]
            y = self.mesh.face_centers[face_i][1]
            z = self.mesh.face_centers[face_i][2]

            values[i] = self.T(x, y, z, t)

        return values

class HeatSource(cedar.Source):
    def __init__(self, Qdot: float | dict[str, float] | Callable = None, volQdot: float | dict[str, float] | Callable = None, shape_func: Callable = None):
        self.Qdot = Qdot
        self.volQdot = volQdot
        self.shape_func = shape_func

        self.mesh: cedar.Mesh3D = None

        self.Qdot_func: Callable = None
        self._value: np.ndarray = None
        self.w: np.ndarray = None

    def initialize(self, t_start: float):
        self.w = self._get_weighting()

        # NOW CALCULATE VALUE USING SHAPE AND QDOT
        self._value = np.zeros(self.mesh.N_cells, dtype = np.float64)

        # If volQdot is given, convert it to Qdot
        if self.volQdot is not None:
            if isinstance(self.volQdot, dict):
                self.Qdot = {}

                for region in self.mesh.regions:
                    self.Qdot[region] = self.volQdot[region] * self.mesh.region_vol[region]

            elif callable(self.volQdot):
                self.Qdot_func = lambda t: sum(self.mesh.region_vol.values())*self.volQdot(t)

            else:
                self.Qdot = self.volQdot * sum(self.mesh.region_vol.values())

        if isinstance(self.Qdot, dict):
            for region in self.mesh.regions:
                if region not in self.Qdot[region]:
                    self.Qdot[region] = 0

                region_i = self.mesh.region_i[region]

                w = self.w[region_i]
                w = w / np.sum(w)

                self._value[region_i] = self.Qdot[region] * w

        elif callable(self.Qdot):
            self.Qdot_func = self.Qdot

        else:
            self._value = self.Qdot * self.w

        self.step(t_start)

    def iterate(self):
        return np.copy(self._value)

    def step(self, t: float):
        if self.Qdot_func is not None:
            self._value = self.Qdot_func(t)*self.w

    def _get_weighting(self):
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
    def __init__(
            self,
            mesh: cedar.Mesh3D = None,
            materials: dict[str, cedar.Material] = {},
            bc: dict[str, cedar.BC] = {},
            source: HeatSource = HeatSource(Qdot = 0)
    ):
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

    def cp_by_cell(self, T):
        cp = np.empty(self.mesh.N_cells)

        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            cp[mask] = mat.cp(T[mask])

        return cp

    def initialize(self, t_start: float):
        self._area_over_L: np.ndarray = self.mesh.face_areas / np.sum(self.mesh.face_d, axis=1)

        for boundary in self.mesh.boundaries:
            if boundary not in self.bc:
                self.bc[boundary] = AdiabaticBC()

        self._rho = self.rho_by_cell()
        self._mass = self._rho*self.mesh.cell_vols

    def iterate(self, res_reduc: float, dt: float) -> float:
        T = self.T.as_continuous_cell_value()

        Qdot = self.source.iterate()
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

    def k_by_cell(self, T):
        k = np.empty(self.mesh.N_cells)
        
        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            k[mask] = mat.k(T[mask])

        return k

    def k_by_face(self, k_cell):
        k_face = np.empty(self.mesh.N_faces)

        f = self.mesh.face_is_interior
        c1 = self.mesh.face_cells_i[f, 0]
        c2 = self.mesh.face_cells_i[f, 1]
        d1 = self.mesh.face_d[f, 0]
        d2 = self.mesh.face_d[f, 1]

        k_face[f] = k_cell[c1]*k_cell[c2] / ((k_cell[c1]*d2 + k_cell[c2]*d1)/(d1+d2))

        # boundary faces
        bf = ~f
        k_face[bf] = k_cell[self.mesh.face_cells_i[bf, 0]]

        return k_face

    def rho_by_cell(self):
        rho = np.empty(self.mesh.N_cells)
        
        for region, mask in self.mesh.region_i.items():
            mat = self.materials[region]
            rho[mask] = mat.rho_rt()

        return rho

    def _A_b(self, k_face, Qdot, cp, dt):
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

            LHS, RHS = self.bc[boundary].iterate(T_boundary, k_boundary)

            rows.extend(c1); cols.extend(c1); data.extend(LHS)
            b[c1] += RHS

        A = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

        return A, b
