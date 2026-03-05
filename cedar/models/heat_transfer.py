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
        Initialize the boundary condition.
        """
        self.boundary: str = None
        self.mesh: cedar.Mesh = None

    def _interpret_user_bc(self, var: np.ndarray | Callable | cedar.Transfer):
        """
        Convert user bc input to a Callable f(x, y, z, t)

        :param var: User input float, np.array, Callable, or cedar.Transfer
        """

        if callable(var):
            return var
        
        elif isinstance(var, cedar.Transfer):
            return lambda x, y, z, t: var.get()
        
        # Must be a number or array, convert to numpy
        var = np.array(var, dtype = np.float64)

        if var.size == 1:
            return lambda x, y, z, t: np.full(self.mesh.boundary_N[self.boundary], var, np.float64)

        else:
            return lambda x, y, z, t: var

    def evaluate_func(self, f: Callable, t: float) -> np.ndarray:
        """
        Evaluate a func f(x, y, z, t) at all face centers on boundary.

        :param f: Function to evaluate.
        :param t: Current simulation time.
        """
        face_i = self.mesh.boundary_i[self.boundary]

        x = self.mesh.face_centers[face_i, 0]
        y = self.mesh.face_centers[face_i, 1]
        z = self.mesh.face_centers[face_i, 2]

        return np.array(f(x, y, z, t), dtype = np.float64)

    @abstractmethod
    def boundary_J(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Compute the heat flux through boundary faces.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Heat flux through boundary faces.
        """
        pass

    @abstractmethod
    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        """
        Compute the temperature at boundary faces.

        :param T_cell: Adjacent cell temperatures.
        :param k: Face thermal conductivity.
        :return: Boundary temperature values.
        """
        pass

    @abstractmethod
    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Return contributions for A, b matrices.

        :param T: Boundary temperature values.
        :param k: Face thermal conductivity.
        :return: (LHS, RHS) diagonal and source term contributions.
        """
        pass

class AdiabaticBC(HeatTransferBC):
    """
    Adiabatic (zero-flux) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces zero normal heat flux across the boundary.
    """
    def __init__(self):
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

    def boundary_J(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return np.zeros_like(T_cell)

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return np.copy(T_cell)

    def initialize(self, t_start: float):
        pass

    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        LHS = np.zeros(self.mesh.boundary_N[self.boundary], dtype = np.float64)
        RHS = np.zeros(self.mesh.boundary_N[self.boundary], dtype = np.float64)

        return LHS, RHS

    def update(self, t: float):
        pass

class FixedTemperatureBC(HeatTransferBC):
    """
    Fixed-temperature (Dirichlet) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces a prescribed boundary temperature, optionally varying in space
    and time.
    """
    def __init__(self, T: float | np.ndarray | Callable):
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

        self.T: float | np.ndarray | Callable = T
        self._T_func: Callable = None

        self.area_over_d: np.ndarray = None

    def boundary_J(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        face_i = self.mesh.boundary_i[self.boundary]
        return k*(T_cell-self.T)/self.mesh.face_d[face_i, 0]

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return self.T

    def initialize(self, t_start: float):
        self._T_func = self._interpret_user_bc(self.T)
        
        face_i = self.mesh.boundary_i[self.boundary]
        self.area_over_d = self.mesh.face_areas[face_i] / self.mesh.face_d[face_i, 0]

    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        LHS = k*self.area_over_d
        RHS = k*self.area_over_d * self.T

        return LHS, RHS

    def update(self, t: float):
        self.T = self.evaluate_func(self._T_func, t)
    
class FixedFluxBC(HeatTransferBC):
    """
    Fixed-flux (Neumann) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces a prescribed boundary flux, optionally varying in space
    and time.
    """
    def __init__(self, J: float | Callable):
        self.boundary: str = None
        self.mesh: cedar.Mesh3D = None

        self.J: float | Callable = J
        self._J_func: Callable = None

        self.area_over_d: np.ndarray = None

    def boundary_J(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        return self.J

    def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
        face_i = self.mesh.boundary_i[self.boundary]
        return T_cell - self.mesh.face_d[face_i, 0]*self.J/k

    def initialize(self, t_start: float):
        self._J_func = self._interpret_user_bc(self.J)

    def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        face_i = self.mesh.boundary_i[self.boundary]
        LHS = np.zeros_like(T)
        RHS = -self.J*self.mesh.face_areas[face_i]

        return LHS, RHS

    def update(self, t: float):
        self.J = self.evaluate_func(self._J_func, t)

# class ConvectiveBC(HeatTransferBC):
#     def __init__(self, T_flow: float | Callable | cedar.FieldView, h: float | Callable | cedar.FieldView):
#         self.boundary: str = None
#         self.mesh: cedar.Mesh3D = None

#         self.T_flow: float | Callable | cedar.FieldView = T_flow
#         self._T_flow_func: Callable = None

#         self.h: float | Callable | cedar.FieldView = h
#         self._h_func: Callable = None

#     def boundary_J(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
#         face_i = self.mesh.boundary_i[self.boundary]
#         d = self.mesh.face_d[face_i, 0]

#         return k*self.h*(T_cell - self.T_flow)/(d*self.h+k)

#     def boundary_T(self, T_cell: np.ndarray, k: np.ndarray) -> np.ndarray:
#         face_i = self.mesh.boundary_i[self.boundary]
#         d = self.mesh.face_d[face_i, 0]
#         return (k*T_cell + self.h*d*self.T_flow)/(self.h*d+k)

#     def initialize(self, t_start: float):
#         if isinstance(self.T_flow, cedar.FieldView):
#             N = self.mesh.boundary_N[self.boundary]
#             self.T_flow.N

#             w = np.zeros((N, self.T_flow.N))

#             for i in range(N):

            

#         self._T_flow_func = self.convert_var_to_func(self.T_flow)
#         self._h_func = self.convert_var_to_func(self._h_func)

#     def matrix_contribution(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
#         face_i = self.mesh.boundary_i[self.boundary]
#         LHS = np.zeros_like(T)
#         RHS = -self.mesh.face_areas[face_i]*self.h*(T-self.T_flow)

#         return LHS, RHS

#     def update(self, t: float):
#         self.T_flow = self.evaluate_func(self._T_flow_func, t)
#         self.h = self.evaluate_func(self._h_func, t)

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
            shape: np.ndarray | Callable = None
        ):
        """
        Initialize the heat source.

        :param Qdot: Total heat input (global, per-region, or time-dependent).
        :param volQdot: Volumetric heat input (global, per-region, or time-dependent).
        :param shape: Spatial weighting function f(x, y, z) or per-region mapping.
        """
        self.Qdot: float | dict[str, float] | Callable | dict[str, Callable] | cedar.FieldView = Qdot
        self._Qdot_func: dict[str, Callable] = None

        self.volQdot = volQdot

        self.shape: Callable | dict[str, Callable] | cedar.FieldView = shape
        self._shape_func = None

        self.mesh: cedar.Mesh3D = None

    def get(self):
        return np.copy(self.Qdot)

    def initialize(self, t_start: float):
        """
        Initialize the source distribution and evaluate initial value.

        :param t_start: Simulation start time.
        """
        self._Qdot_func = {}
        self._shape_func = {}

        if self.volQdot is not None:
            self.Qdot = self._interpret_user_volQdot(self.volQdot)

        self._Qdot_func = self._interpret_user_Qdot(self.Qdot)
        self._shape_func = self._interpret_user_shape(self.shape)

        self.Qdot = np.zeros(self.mesh.N_cells, dtype = np.float64)

    def update(self, t: float):
        for region in self.mesh.regions:
            region_i = self.mesh.region_i[region]

            self.Qdot[region_i] = self._Qdot_func[region](t)*self._shape_func[region]()

    def _interpret_user_Qdot(self, user_Qdot) -> dict[str, Callable]:
        Qdot = {}

        if isinstance(user_Qdot, dict):
            for region in self.mesh.regions:
                if region not in user_Qdot:
                    user_Qdot[region] = 0

        else:
            global_Qdot = user_Qdot
            user_Qdot = {}
            total_vol = sum(self.mesh.region_vol.values())
            total_vol = np.array(total_vol, dtype = np.float64)

            for region in self.mesh.regions:
                region_fraction = self.mesh.region_vol[region] / total_vol

                if callable(global_Qdot):
                    user_Qdot[region] = lambda t: global_Qdot(t) * region_fraction

                else:
                    user_Qdot[region] = global_Qdot * region_fraction

        for region in user_Qdot:
            region_Qdot = user_Qdot[region]

            if callable(region_Qdot):
                Qdot[region] = region_Qdot
            
            else:
                Qdot[region] = lambda t: region_Qdot

        return Qdot

    def _interpret_user_shape(self, user_shape) -> dict[str, Callable]:
        shape = {}

        if isinstance(user_shape, dict):
            for region in self.mesh.regions:
                if region not in user_shape:
                    user_shape[region] = None

        else:
            global_shape = user_shape

            user_shape = {}
            for region in self.mesh.regions:
                user_shape[region] = global_shape

        for region in user_shape:
            region_i = self.mesh.region_i[region]

            region_shape = user_shape[region]

            if callable(region_shape):
                x = self.mesh.cell_centers[region_i, 0]
                y = self.mesh.cell_centers[region_i, 1]
                z = self.mesh.cell_centers[region_i, 2]

                values = region_shape(x, y, z)
                
            else:
                values = np.ones(self.mesh.region_N[region], dtype = np.float64)
            
            values *= self.mesh.cell_vols[region_i]
            values = values / sum(values)

            shape[region] = lambda: values

        return shape

    def _interpret_user_volQdot(self, volQdot) -> dict:
        Qdot = {}

        if isinstance(volQdot, dict):
            for region in self.mesh.regions:
                if region not in volQdot:
                    volQdot[region] = 0

        else:
            global_volQdot = volQdot

            volQdot = {}

            for region in self.mesh.regions:
                volQdot[region] = global_volQdot
    
        for region in volQdot:
            if callable(volQdot[region]):
                Qdot[region] = lambda t: self.mesh.region_vol[region] * volQdot[region](t)

            else:
                Qdot[region] = self.volQdot[region] * self.mesh.region_vol[region]

        return Qdot

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
        :param source: Heat generation source.
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

        self.J: cedar.Field = self.add_field("J", is_on_regions = False, is_on_boundaries = True)

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
        self._area_over_d: np.ndarray = self.mesh.face_areas / np.sum(self.mesh.face_d, axis=1)

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
        T = self.T.get()
        T_ic = self.T.get_ic()

        Qdot = self.source.get()
        volQdot = Qdot / self.mesh.cell_vols

        k_cell = self.k_by_cell(T)
        k_face = self.k_by_face(k_cell)
        cp = self.cp_by_cell(T)

        A, b = self._A_b(k_face, Qdot, cp, T_ic, dt)

        res = cedar.functions.residual(A, b, T)
        
        tgt_res = res_reduc * res

        T, _ = cg(A, b, T, rtol = tgt_res)

        self.T.set(T)
        self.Qdot.set(Qdot)
        self.volQdot.set(volQdot)
        self.k.set(k_cell)
        self.rho.set(self._rho)
        self.cp.set(cp)

        for boundary in self.mesh.boundaries:
            face_i = self.mesh.boundary_i[boundary]
            T_cell = T[self.mesh.face_cells_i[face_i, 0]]

            T_boundary_face = self.bc[boundary].boundary_T(T_cell, k_face[face_i])
            self.T.set(T_boundary_face, boundary)

            J_boundary_face = self.bc[boundary].boundary_J(T_cell, k_face[face_i])
            self.J.set(J_boundary_face, boundary)

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

    def _A_b(self, k_face, Qdot, cp, T_ic, dt):
        """
        Assemble the linear system A T = b.

        :param k_face: Face thermal conductivity.
        :param Qdot: Cell heat source.
        :param cp: Cell specific heat capacity.
        :param T_ic: Value of T at initial condition.
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

            rows.extend(range(N)); cols.extend(range(N)); data.extend(C)

            b += T_ic * C

        # -----------------------------------------
        #  INTERIOR FACES
        # -----------------------------------------
        face_data = zip(
            self.mesh.face_cells_i[self.mesh.face_is_interior],
            k_face[self.mesh.face_is_interior],
            self._area_over_d[self.mesh.face_is_interior]
        )

        for cell_ids, k, area_over_d in face_data:
            C = k * area_over_d
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