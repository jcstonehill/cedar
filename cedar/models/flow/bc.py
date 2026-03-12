from typing import Callable

import cedar


class InletBC(cedar.Term):
    """
    Inlet boundary condition for :class:`cedar.Flow`.

    Enforces a prescribed stagnation temperature and mass flow rate.

    Parameters
    ----------
    T0 : float or f(t), optional
        Known stagnation temperature to apply to boundary in [K].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. T0 = 300
        - A function f(t) evaluated at simulation time. ex. T0 = lambda t: t + 300

    mdot : float or f(t), optional
        Known mass flow rate to apply to boundary in [kg/s].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. mdot = 1
        - A function f(t) evaluated at simulation time. ex. mdot = lambda t: t + 1

    Attributes
    ----------
    T0 : float or f(t)
        Known stagnation temperature to apply to boundary in [K].

    mdot : float or f(t)
        Known mass flow rate to apply to boundary in [kg/s].

    """

    def __init__(
        self,
        T0: float | Callable | cedar.Field = None,
        mdot: float | Callable | cedar.Field = None,
    ):
        self.T0: float | Callable | cedar.Field = T0
        self.mdot: float | Callable | cedar.Field = mdot

    def _initialize(self):
        T0 = self.T0
        mdot = self.mdot

        if callable(T0):
            T0 = lambda x, y, z, t: T0(t)

        if callable(mdot):
            mdot = lambda x, y, z, t: mdot(t)

        self._add_value("T0", T0)
        self._add_value("mdot", mdot)


class OutletBC(cedar.Term):
    """
    Outlet boundary condition for :class:`cedar.Flow`.

    Enforces a prescribed stagnation pressure.

    Parameters
    ----------
    P0 : float or f(t), optional
        Known stagnation pressure to apply to boundary in [Pa].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. P0 = 101325
        - A function f(t) evaluated at simulation time. ex. P0 = lambda t: t + 101325

    Attributes
    ----------
    P0 : float or f(t)
        Known stagnation pressure to apply to boundary in [Pa].

    """

    def __init__(self, P0: float | Callable | cedar.Field = None):
        self.P0: float | Callable | cedar.Field = P0

    def _initialize(self):
        P0 = self.P0

        if callable(P0):
            P0 = lambda x, y, z, t: P0(t)

        self._add_value("P0", P0)
