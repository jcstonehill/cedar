import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve, cg, bicgstab, gmres, cgs
from scipy.sparse import lil_matrix

import cedar


@dataclass
class Params:
    fluid: cedar.base.Fluid = None
    Dh: float = None
    A: float = None
    P_wall: float = None
    eps: float = None

@dataclass
class Vars:
    inlet: cedar.FlowStateVar
    outlet: cedar.FlowStateVar

    Qdot: cedar.ScalarVar
    Qdot_shape: cedar.Mesh1DVar

    T: cedar.Mesh1DVar
    P: cedar.Mesh1DVar
    u: cedar.Mesh1DVar

    rho: cedar.Mesh1DVar
    mu: cedar.Mesh1DVar
    cp: cedar.Mesh1DVar
    k: cedar.Mesh1DVar

    e: cedar.Mesh1DVar
    E: cedar.Mesh1DVar

    Re: cedar.Mesh1DVar
    Pr: cedar.Mesh1DVar

    ff: cedar.Mesh1DVar
    htc: cedar.Mesh1DVar
    Nu: cedar.Mesh1DVar

    T_wall: cedar.Mesh1DVar

class Flow(cedar.base.Model):
    """
    Solves the Compressible Euler Eqs.

    Used to calculate bulk flow temperature and be coupled with a ``Thermal``
    model to calculate heat removal.
    """
    @property
    def mesh(self) -> cedar.Mesh1D:
        return self._mesh
    
    @property
    def params(self) -> Params:
        return self._params
    
    @property
    def vars(self) -> Vars:
        return self._vars

    def __init__(self, name, mesh: cedar.Mesh1D, fluid: cedar.base.Fluid,
                 Dh: float, A: float = None, P_wall: float = None, eps: float = 0):
        self.name = name
        self._mesh = mesh

        A = 0.25*np.pi*Dh**2 if A is None else A
        P_wall = np.pi*Dh if P_wall is None else P_wall

        self._params = Params(fluid, Dh, A, P_wall, eps)

        self._vars = Vars(
            inlet = self.add_var(cedar.FlowStateVar(self.name, "inlet")),
            outlet = self.add_var(cedar.FlowStateVar(self.name, "outlet")),

            Qdot = self.add_var(cedar.ScalarVar(self.name, "Q_dot", 0, "W", self.mesh)),
            Qdot_shape = self.add_var(cedar.Mesh1DVar(self.name, "Qdot_shape", 1, None, self.mesh)),
            
            T = self.add_var(cedar.Mesh1DVar(self.name, "T", 300, "K", self.mesh)),
            P = self.add_var(cedar.Mesh1DVar(self.name, "P", 101325, "Pa", self.mesh)),
            u = self.add_var(cedar.Mesh1DVar(self.name, "u", 1, "m/s", self.mesh)),
            rho = self.add_var(cedar.Mesh1DVar(self.name, "rho", 1, "kg/m3", self.mesh)),
            mu = self.add_var(cedar.Mesh1DVar(self.name, "mu", 1, "Pa-s", self.mesh)),
            cp = self.add_var(cedar.Mesh1DVar(self.name, "cp", 1, "J/kg-K", self.mesh)),
            k = self.add_var(cedar.Mesh1DVar(self.name, "k", 1, "W/m-K", self.mesh)),
            
            e = self.add_var(cedar.Mesh1DVar(self.name, "e", 1, "J/kg", self.mesh)),
            E = self.add_var(cedar.Mesh1DVar(self.name, "E", 1, "J/kg", self.mesh)),

            Re = self.add_var(cedar.Mesh1DVar(self.name, "Re", 1, None, self.mesh)),
            Pr = self.add_var(cedar.Mesh1DVar(self.name, "Pr", 1, None, self.mesh)),

            ff = self.add_var(cedar.Mesh1DVar(self.name, "ff", 1, None, self.mesh)),
            htc = self.add_var(cedar.Mesh1DVar(self.name, "htc", 1, None, self.mesh)),
            Nu = self.add_var(cedar.Mesh1DVar(self.name, "Nu", 1, None, self.mesh)),

            T_wall = self.add_var(cedar.Mesh1DVar(self.name, "T_wall", 300, "K", self.mesh))
        )

    def check(self):
        pass

    def setup(self):
        self.vars.T.set(self.vars.T.initial)
        self.vars.P.set(self.vars.P.initial)
        self.vars.u.set(self.vars.u.initial)

    def iterate(self, t0 = None, dt = None, res_reduc = 1e-2):
        # Get some params
        A = self.params.A
        P_wall = self.params.P_wall
        eps = self.params.eps

        Qdot = self.vars.Qdot.val
        Qdot_shape = self.vars.Qdot_shape.val

        if sum(Qdot_shape) == 0:
            Qdot_shape = np.ones(self.mesh.N_cells)

        Qdot_shape = np.insert(Qdot_shape, 0, 0)
        
        Qdot_shape = Qdot_shape / sum(Qdot_shape)
        Qdot = Qdot*Qdot_shape

        # Initialize output arrays
        T = cedar.helper.add_upwind_ghost_cell(self.vars.T.val)
        P = cedar.helper.add_upwind_ghost_cell(self.vars.P.val)
        u = cedar.helper.add_upwind_ghost_cell(self.vars.u.val)
        E = cedar.helper.add_upwind_ghost_cell(self.vars.E.val)
        
        rho, k, cp, mu = self._properties(T, P)

        for i in range(100):
            rho_prev = np.copy(rho)
            u_prev = np.copy(u)
            T_prev = np.copy(T)
            P_prev = np.copy(P)

            Re = self._Re(rho, u, self.params.Dh, mu)

            Pr = self._Pr(k, cp, mu)
            ff = self._ff(eps, self.params.Dh, Re)

            u, mass_res = self._cons_mass(u, rho, A, dt)
            P, mom_res = self._cons_momentum(P, T, rho, u, ff, self.params.Dh, dt)
            E, energy_res = self._cons_energy(E, T, P, rho, u, Qdot, A, dt)

            e = self._e(E, u)
            T = self._T(e, P, u)

            rho, k, cp, mu = self._properties(T, P)

            # res = np.linalg.norm(np.concatenate((rho-rho_prev, u-u_prev, T-T_prev, P-P_prev)))
            res = np.max([mass_res, mom_res, energy_res])

            if res <= 1e-6:
                self._post(T, P, u, rho, k, cp, mu, e, E, Qdot, Re, Pr, ff, A)
                return res
            
            

        self.log_error("Flow model did not converge.")

            # rho, k, cp, mu = self._properties(T, P)
            # #print("Re", Re, "Pr", Pr, "ff", ff, "u", u, "P", P, "E", E, "e", e, "T", T, "rho", rho)
        
            # #`print`(res)
            # if res <= 1e-6:
            #     self._post(T, P, u, rho, k, cp, mu, e, E, Qdot, Re, Pr, ff, A)
                
            #     return res

    def _ff(self, eps, Dh, Re):
        return cedar.helper.churchill(eps, Dh, Re)

    def _cons_energy(self, E, T, P, rho, u, Qdot, area, dt):
        N = self.mesh.N_cells + 1

        # Sparse matrix assembly (LIL is efficient for incremental filling)
        A = lil_matrix((N, N))
        b = np.zeros(N)

        e_in = self.params.fluid.e_from_T_P(T[0], P[0])
        E_in = e_in + 0.5 * u[0]**2

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

        # if dt is not None:
        #     for i in range(1, self.mesh.N_cells+2):
        #         A[i][i] += rho[i]*self.mesh.dz/dt
        #         b[i] += self.vars.rho.initial[i]*self.vars.E.initial[i]*self.mesh.dz/dt

        # Convert once to CSR for fast solving
        A = A.tocsr()

        res = cedar.helper.residual(A, b, E)

        E = spsolve(A, b)

        return E, res

    def _cons_mass(self, u, rho, area, dt):
        N = self.mesh.N_cells + 1

        A = lil_matrix((N, N))
        b = np.zeros(N)

        # Boundary condition
        A[0, 0] = 1.0
        b[0] = self.vars.inlet.val[2] / (rho[0] * area)

        # Interior equations
        for i in range(1, N):
            A[i, i-1] =  rho[i-1]
            A[i, i]   = -rho[i]

        # if dt is not None:
        #     for i in range(1, self.mesh.N_cells+2):
        #         b[i] += (rho[i] - self.vars.rho.initial[i])*self.mesh.dz/dt

        # Convert to CSR for efficient solve
        A = A.tocsr()

        res = cedar.helper.residual(A, b, u)

        u = spsolve(A, b)

        return u, res

    def _cons_momentum(self, P, T, rho, u, ff, Dh, dt):
        N = self.mesh.N_cells + 1

        A = lil_matrix((N, N))
        b = np.zeros(N)

        # Outlet boundary condition
        A[-1, -1] = 1.0
        b[-1] = cedar.helper.P0_to_P(
            self.vars.outlet.val[1],
            T[-1],
            u[-1],
            self.params.fluid
        )

        # Interior equations (backward sweep)
        for i in range(self.mesh.N_cells - 1, -1, -1):
            A[i, i]   = 1.0
            A[i, i+1] = -1.0

            dP_acc = rho[i+1] * u[i+1]**2 - rho[i] * u[i]**2
            dP_fric = (ff[i] * self.mesh.dz * rho[i] * u[i]**2) / (2 * Dh)

            b[i] = dP_acc + dP_fric

        # Convert once to CSR for efficient solve
        A = A.tocsr()

        # if dt is not None:
            #     for i in range(self.mesh.N_cells):
            #         b[i] += (rho[i]*u[i] - self.vars.rho.initial[i]*self.vars.u.initial[i])*self.mesh.dz/dt

        res = cedar.helper.residual(A, b, P)

        P = spsolve(A, b)

        return P, res

    def _e(self, E, u):
        return E - 0.5*u**2

    def _htc(self, Nu, k, d):
        return Nu*k/d

    def _Nu(self, Re, Pr, T_wall, T, z, d):
        return cedar.helper.westinghouse_modified_mccarthy_wolf(Re, Pr, T_wall, T, z, d)
    
    def _post(self, T, P, u, rho, k, cp, mu, e, E, Qdot, Re, Pr, ff, A):
        T0_in = self.vars.inlet.val[0]
        P0_in = cedar.helper.P_to_P0(P[0], u[0], rho[0])
        mdot_in = self.vars.inlet.val[2]

        T0_out = cedar.helper.T_to_T0(T[-1], u[-1], cp[-1])
        P0_out = self.vars.outlet.val[1]
        mdot_out = rho[-1]*u[-1]*A

        z = np.array([pos[2] for pos in self.mesh.cell_centers])
        T_wall, Nu, htc = self._wall_heating_vals(Qdot[1:], Re[1:], Pr[1:], T[1:], k[1:], z, self.params.P_wall, self.params.Dh)

        self.vars.inlet.set([T0_in, P0_in, mdot_in])
        self.vars.outlet.set([T0_out, P0_out, mdot_out])

        self.vars.T.set(T[1:])
        self.vars.P.set(P[1:])
        self.vars.u.set(u[1:])

        self.vars.rho.set(rho[1:])
        self.vars.mu.set(mu[1:])
        self.vars.cp.set(cp[1:])
        self.vars.k.set(k[1:])

        self.vars.e.set(e[1:])
        self.vars.E.set(E[1:])

        self.vars.Re.set(Re[1:])
        self.vars.Pr.set(Pr[1:])

        self.vars.ff.set(ff[1:])
        self.vars.htc.set(htc)
        self.vars.Nu.set(Nu)

        self.vars.T_wall.set(T_wall)

    def _Pr(self, k, cp, mu):
        return cedar.helper.prandtl(k, cp, mu)

    def _properties(self, T, P):
        return (self.params.fluid.rho_from_T_P(T, P),
                self.params.fluid.k_from_T_P(T, P),
                self.params.fluid.cp_from_T_P(T, P),
                self.params.fluid.mu_from_T_P(T, P))

    def _Re(self, rho, u, Dh, mu):
        return cedar.helper.reynolds(rho, u, Dh, mu)
    
    def _T(self, e, P, u):
        T = self.params.fluid.T_from_e_P(e, P)
        T[0] = cedar.helper.T0_to_T(self.vars.inlet.val[0], P[0], u[0], self.params.fluid)

        return T
    
    def _wall_heating_vals(self, Qdot, Re, Pr, T, k, z, P_wall, Dh):
        A_wall = P_wall*self.mesh.dz

        T_wall = T

        for _ in range(1000):
            T_wall_prev = np.copy(T_wall)
            Nu = self._Nu(Re, Pr, T_wall, T, z, Dh)
            htc = self._htc(Nu, k, Dh)

            T_wall = Qdot/(htc*A_wall) + T

            if np.linalg.norm(T_wall-T_wall_prev) <= 1:
                return T_wall, Nu, htc
            
        self.log_error("Walling heating values did not converge.")