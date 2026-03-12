import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix

import cedar


class Flow(cedar.Model):
    """
    Flow model using FVM.

    Solves the compressible euler flow equations on a 1D mesh.

    Parameters
    ----------
    name : str
        Name of this flow channel. The created Mesh1D will have 1 region with this
        name and boundaries named "{name}:in" and "{name}:out" for the inlet and
        outlet respectively.
    N_cells : int
        Number of cells in the created Mesh1D.
    L : float
        Length of flow channel in [m].
    start : tuple(float, float, float), optional
        Position of inlet in (x, y, z) coordinates in [m].
    basis : {"x", "y", "z", "-x", "-y", "-z"}, optional
        The direction of the flow. Defaults to (0, 0, 0).
    A : float, optional
        Cross sectional flow area in [m^2].
    P_wall : float, optional
        Perimeter of wall in [m].
    eps : float, optional
        Wall roughness in [m].
    fluid : cedar.Fluid, optional
        Fluid properties object.

    Attributes
    ----------
    mesh : cedar.Mesh1D
        The mesh that the flow equations will be solved on.
    A : float
        Cross sectional flow area in [m^2].
    P_wall : float
        Perimeter of wall in [m].
    eps : float
        Wall roughness in [m].
    fluid : cedar.Fluid
        Fluid properties object.
    bcs : dict[str, cedar.Term]
        Boundary Conditions. The keys are boundary names in Flow.mesh.
    sources : dict[str, cedar.FlowSource]
        Source Terms. The key is the region name in Flow.mesh.
    inlet : cedar.InletBC
        Inlet boundary condition.
    outlet : cedar.OutletBC
        Outlet boundary condition.
    T : cedar.Field
        Static temperature in [K].
    P : cedar.Field
        Static pressure in [Pa].
    u : cedar.Field
        Flow velocity in [m/s].
    rho : cedar.Field
        Mass density in [kg/m^3].
    mu : cedar.Field
        Dynamic viscosity in [Pa-s].
    cp : cedar.Field
        Specific heat capacity in [J/kg-K].
    k : cedar.Field
        Thermal conductivity in [W/m-K].
    e : cedar.Field
        Specific static internal energy in [J/kg].
    E : cedar.Field
        Specific stagnation internal energy in [J/kg].
    Re : cedar.Field
        Reynolds number.
    Pr : cedar.Field
        Prandtl number.
    ff : cedar.Field
        Friction factor.
    Nu : cedar.Field
        Nusselt number.
    htc : cedar.Field
        Heat transfer coefficient in [W/m^2-K].
    Qdot : cedar.Field
        Heat transferred to the flow in [W].

    """

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
        basis="z",
        A: float = None,
        P_wall: float = None,
        eps: float = None,
        fluid: cedar.Fluid = None,
    ):
        self._mesh: cedar.Mesh1D = cedar.Mesh1D(
            N_cells=N_cells,
            L=L,
            region=name,
            start_boundary=name + ":in",
            end_boundary=name + ":out",
            start=start,
            basis=basis,
        )
        self.A = A
        self.P_wall = P_wall
        self.eps = eps
        self.fluid = fluid

        self.bcs: dict[str, cedar.Term] = {}
        self.sources: dict[str, cedar.FlowSource] = {}

        self.inlet = cedar.InletBC()
        self.outlet = cedar.OutletBC()

        self.T = self._add_field("T", "Flow", 300)
        self.P = self._add_field("P", "Flow", 101325)
        self.u = self._add_field("u", "Flow", 1)

        self.rho = self._add_field("rho", "Flow", 1)
        self.mu = self._add_field("mu", "Flow", 1)
        self.cp = self._add_field("cp", "Flow", 1)
        self.k = self._add_field("k", "Flow", 1)

        self.e = self._add_field("e", "Flow", 1)
        self.E = self._add_field("E", "Flow", 1)

        self.Re = self._add_field("Re", "Flow", 1)
        self.Pr = self._add_field("Pr", "Flow", 1)

        self.ff = self._add_field("ff", "Flow", 0.1)

        self.T_wall = self._add_field("T_wall", "Flow", 300)
        self.Nu = self._add_field("Nu", "Flow", 1)
        self.htc = self._add_field("htc", "Flow", 1)
        self.Qdot = self._add_field("Qdot", "Flow", 0)

    def add_source(self, source: cedar.FlowSource):
        """
        Add source term to flow model.

        Since heat transfer can only come from the wall, only one source term
        can be defined. Adding a second source term will override the previous.

        Parameters
        ----------
        source : cedar.FlowSource
            Source of heat transfer to flow.

        """
        self.sources[self.mesh.regions[0]] = source

    def _build(self):
        if not self.sources:
            self.sources: dict[str, cedar.FlowSource] = {
                self.mesh.regions[0]: cedar.FlowQdotSource(0)
            }

        self.bcs: dict[str, cedar.InletBC | cedar.OutletBC] = {
            self.mesh.boundaries[0]: self.inlet,
            self.mesh.boundaries[1]: self.outlet,
        }

    def _initialize(self):
        self.L = np.linalg.norm(self.mesh.end - self.mesh.start)
        self.Dh = 4 * self.A / self.P_wall

    def _iterate(self, res_reduc: float, dt: float) -> float:
        T = cedar.functions.add_upwind_ghost_cell(self.T.get())
        P = cedar.functions.add_upwind_ghost_cell(self.P.get())
        u = cedar.functions.add_upwind_ghost_cell(self.u.get())
        E = cedar.functions.add_upwind_ghost_cell(self.E.get())

        rho, k, cp, mu = self._properties(T, P)

        for _ in range(100):
            rho_prev = np.copy(rho)

            Re = cedar.functions.reynolds(rho, u, self.L, mu)
            Pr = cedar.functions.prandtl(k, cp, mu)
            ff = cedar.functions.churchill(self.eps, self.Dh, Re)

            Qdot, T_wall, Nu, htc = self.sources[
                self.mesh.regions[0]
            ]._wall_heating_values(
                T[1:], Re[1:], Pr[1:], k[1:], self.Dh, self.P_wall, self.mesh.ds
            )

            Qdot = np.insert(Qdot, 0, 0)

            u, mass_res = self._conserve_mass(u, rho, self.A, dt)
            P, mom_res = self._conserve_momentum(P, rho, u, ff, self.Dh, dt)
            E, energy_res = self._conserve_energy(E, P, rho, u, Qdot, self.A, dt)

            e = E - 0.5 * u**2

            T = self.fluid.T_from_e_P(e, P)

            rho, k, cp, mu = self._properties(T, P)

            rho_res = np.linalg.norm(rho - rho_prev)

            res = np.max([mass_res, mom_res, energy_res, rho_res])

            if res <= 1e-6:
                self.T._set(T[1:])
                self.P._set(P[1:])
                self.u._set(u[1:])

                self.rho._set(rho[1:])
                self.mu._set(mu[1:])
                self.cp._set(cp[1:])
                self.k._set(k[1:])

                self.e._set(e[1:])
                self.E._set(E[1:])

                self.Re._set(Re[1:])
                self.Pr._set(Pr[1:])

                self.ff._set(ff[1:])

                self.T_wall._set(T_wall)
                self.Nu._set(Nu)
                self.htc._set(htc)

                self.Qdot._set(Qdot[1:])

                return res

    def _conserve_energy(self, E, P, rho, u, Qdot, area, dt):
        N = self.mesh.N_cells + 1

        # Sparse matrix assembly (LIL is efficient for incremental filling)
        A = lil_matrix((N, N))
        b = np.zeros(N)

        P0_in = cedar.functions.P_to_P0(P[0], u[0], rho[0])
        T0_in = self.bcs[self.mesh.boundaries[0]].values["T0"][0]
        E_in = self.fluid.e_from_T_P(T0_in, P0_in)

        # Boundary condition
        A[0, 0] = 1.0
        b[0] = E_in

        # Interior equations
        for i in range(1, N):
            A[i, i] = u[i] * rho[i]
            A[i, i - 1] = -u[i - 1] * rho[i - 1]

            b[i] = u[i - 1] * P[i - 1] - u[i] * P[i] + Qdot[i] / area

        A = A.tocsr()
        res = cedar.functions.residual(A, b, E)
        E = spsolve(A, b)

        return E, res

    def _conserve_mass(self, u, rho, area, dt):
        N = self.mesh.N_cells + 1

        A = lil_matrix((N, N))
        b = np.zeros(N)

        A[0, 0] = 1.0
        b[0] = self.bcs[self.mesh.boundaries[0]].values["mdot"][0] / (rho[0] * area)

        # Interior equations
        for i in range(1, N):
            A[i, i - 1] = rho[i - 1]
            A[i, i] = -rho[i]

        A = A.tocsr()
        res = cedar.functions.residual(A, b, u)
        u = spsolve(A, b)

        return u, res

    def _conserve_momentum(self, P, rho, u, ff, Dh, dt):
        N = self.mesh.N_cells + 1

        A = lil_matrix((N, N))
        b = np.zeros(N)

        P0_out = self.bcs[self.mesh.boundaries[1]].values["P0"][0]

        A[-1, -1] = 1.0
        b[-1] = cedar.functions.P0_to_P(P0_out, u[-1], rho[-1])

        for i in range(self.mesh.N_cells - 1, -1, -1):
            A[i, i] = 1.0
            A[i, i + 1] = -1.0

            dP_acc = rho[i + 1] * u[i + 1] ** 2 - rho[i] * u[i] ** 2
            dP_fric = (ff[i] * self.mesh.ds * rho[i] * u[i] ** 2) / (2 * Dh)

            b[i] = dP_acc + dP_fric

        A = A.tocsr()
        res = cedar.functions.residual(A, b, P)
        P = spsolve(A, b)

        return P, res

    def _properties(self, T, P):
        return (
            self.fluid.rho_from_T_P(T, P),
            self.fluid.k_from_T_P(T, P),
            self.fluid.cp_from_T_P(T, P),
            self.fluid.mu_from_T_P(T, P),
        )
