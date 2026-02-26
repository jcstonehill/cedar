import numpy as np
from scipy.sparse.linalg import spsolve, cg, bicgstab, gmres, cgs
from scipy.sparse import lil_matrix

import cedar

class InletBC(cedar.BC):
    def __init__(self, T0: float = None, mdot: float = None):
        self.T0 = T0
        self.mdot = mdot

        self.mesh: cedar.Mesh1D = None

    def initialize(self, t_start: float):
        pass

    def step(self, t: float):
        pass

class OutletBC(cedar.BC):
    def __init__(self, P0: float = None):
        self.P0 = P0

        self.mesh: cedar.Mesh1D = None

    def initialize(self, t_start: float):
        pass

    def step(self, t: float):
        pass

class FlowHeatSource(cedar.Source):
    def __init__(self, Qdot: float):
        self.mesh: cedar.Mesh1D = None
        self.Qdot = Qdot

    def initialize(self, t_start: float):
        self.Qdot = np.ones(self.mesh.N_cells)*self.Qdot/self.mesh.N_cells

    def step(self, t: float):
        pass

class Flow(cedar.Model):
    @property
    def mesh(self) -> cedar.Mesh1D:
        return self._mesh
    
    @mesh.setter
    def mesh(self, val):
        cedar.Log.error("Unable to set Flow.mesh after instantiation.")

    def __init__(
            self,
            name: str,
            N_cells: int,
            L: float,
            start: np.ndarray = (0, 0, 0),
            basis = "z",
            A: float = None,
            P_wall: float = None,
            eps: float = None,
            fluid: cedar.Fluid = None,
            source: FlowHeatSource = None
    ):
        self._mesh: cedar.Mesh1D = cedar.Mesh1D(
            N_cells = N_cells,
            L = L, 
            region = name, 
            start_boundary = name+":in",
            end_boundary = name+":out",
            start = start,
            basis = basis
        )
        self.A = A
        self.P_wall = P_wall
        self.eps = eps
        self.fluid = fluid
        self.source: FlowHeatSource = FlowHeatSource(0) if source is None else source

        self.inlet = InletBC()
        self.outlet = OutletBC()

        self.T = self.add_field("T", requires_ic = True)
        self.P = self.add_field("P", requires_ic = True)
        self.u = self.add_field("u", requires_ic = True)

        self.rho = self.add_field("rho")
        self.mu = self.add_field("mu")
        self.cp = self.add_field("cp")
        self.k = self.add_field("k")

        self.e = self.add_field("e")
        self.E = self.add_field("E")

        self.Re = self.add_field("Re")
        self.Pr = self.add_field("Pr")

        self.ff = self.add_field("ff")
        self.htc = self.add_field("htc")
        self.Nu = self.add_field("Nu")

        self.T_wall = self.add_field("T_wall")

    def initialize(self, t_start: float):
        self.Dh = 4*self.A/self.P_wall
        self.L = np.linalg.norm(self.mesh.end-self.mesh.start)

        if self.source is None:
            self.source = FlowHeatSource(0)

        self.bc = {
            self.mesh.boundaries[0] : self.inlet,
            self.mesh.boundaries[1] : self.outlet
        }

    def iterate(self, res_reduc: float, dt: float) -> float:
        T = cedar.functions.add_upwind_ghost_cell(self.T.as_continuous_cell_value())
        P = cedar.functions.add_upwind_ghost_cell(self.P.as_continuous_cell_value())
        u = cedar.functions.add_upwind_ghost_cell(self.u.as_continuous_cell_value())
        E = cedar.functions.add_upwind_ghost_cell(self.E.as_continuous_cell_value())

        # Need to add ghost cell to Qdot as well.
        Qdot = np.insert(self.source.Qdot, 0, 0)

        rho, k, cp, mu = self._properties(T, P)

        for _ in range(100):
            rho_prev = np.copy(rho)

            Re = cedar.functions.reynolds(rho, u, self.L, mu)
            Pr = cedar.functions.prandtl(k, cp, mu)
            ff = cedar.functions.churchill(self.eps, self.Dh, Re)

            u, mass_res = self._conserve_mass(u, rho, self.A, dt)
            P, mom_res = self._conserve_momentum(P, rho, u, ff, self.Dh, dt)
            E, energy_res = self._conserve_energy(E, P, rho, u, Qdot, self.A, dt)

            e = E - 0.5*u**2
            T = self.fluid.T_from_e_P(e, P)

            rho, k, cp, mu = self._properties(T, P)

            rho_res = np.linalg.norm(rho-rho_prev)

            res = np.max([mass_res, mom_res, energy_res, rho_res])

            if res <= 1e-6:
                region = self.mesh.regions[0]
                self.T.values[region] = T[1:]
                self.P.values[region] = P[1:]
                self.u.values[region] = u[1:]
                
                self.rho.values[region] = rho[1:]
                self.mu.values[region] = mu[1:]
                self.cp.values[region] = cp[1:]
                self.k.values[region] = k[1:]
                
                self.e.values[region] = e[1:]
                self.E.values[region] = E[1:]
                
                self.Re.values[region] = Re[1:]
                self.Pr.values[region] = Pr[1:]

                self.ff.values[region] = ff[1:]

                return res

    def _conserve_energy(self, E, P, rho, u, Qdot, area, dt):
        N = self.mesh.N_cells + 1

        # Sparse matrix assembly (LIL is efficient for incremental filling)
        A = lil_matrix((N, N))
        b = np.zeros(N)

        P0_in = cedar.functions.P_to_P0(P[0], u[0], rho[0])
        T0_in = self.bc[self.mesh.boundaries[0]].T0
        E_in = self.fluid.e_from_T_P(T0_in, P0_in)

        # Boundary condition
        A[0, 0] = 1.0
        b[0] = E_in

        # Interior equations
        for i in range(1, N):
            A[i, i]   =  u[i]   * rho[i]
            A[i, i-1] = -u[i-1] * rho[i-1]

            b[i] = (
                u[i-1] * P[i-1]
                - u[i] * P[i]
                + Qdot[i] / area
            )

        A = A.tocsr()
        res = cedar.functions.residual(A, b, E)
        E = spsolve(A, b)

        return E, res

    def _conserve_mass(self, u, rho, area, dt):
        N = self.mesh.N_cells + 1

        A = lil_matrix((N, N))
        b = np.zeros(N)

        A[0, 0] = 1.0
        b[0] = self.bc[self.mesh.boundaries[0]].mdot / (rho[0] * area)

        # Interior equations
        for i in range(1, N):
            A[i, i-1] =  rho[i-1]
            A[i, i]   = -rho[i]

        A = A.tocsr()
        res = cedar.functions.residual(A, b, u)
        u = spsolve(A, b)

        return u, res
    
    def _conserve_momentum(self, P, rho, u, ff, Dh, dt):
        N = self.mesh.N_cells + 1

        A = lil_matrix((N, N))
        b = np.zeros(N)

        P0_out = self.bc[self.mesh.boundaries[1]].P0

        A[-1, -1] = 1.0
        b[-1] = cedar.functions.P0_to_P(P0_out, u[-1], rho[-1])

        for i in range(self.mesh.N_cells - 1, -1, -1):
            A[i, i]   = 1.0
            A[i, i+1] = -1.0

            dP_acc = rho[i+1] * u[i+1]**2 - rho[i] * u[i]**2
            dP_fric = (ff[i] * self.mesh.ds * rho[i] * u[i]**2) / (2 * Dh)

            b[i] = dP_acc + dP_fric

        A = A.tocsr()
        res = cedar.functions.residual(A, b, P)
        P = spsolve(A, b)

        return P, res

    def _properties(self, T, P):
        return (self.fluid.rho_from_T_P(T, P),
                self.fluid.k_from_T_P(T, P),
                self.fluid.cp_from_T_P(T, P),
                self.fluid.mu_from_T_P(T, P))