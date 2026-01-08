import numpy as np

import cedar


def add_upwind_ghost_cell(array: np.ndarray) -> np.ndarray:
    """
    Prepend an upwind ghost cell to a 1D array.

    The ghost cell value is set equal to the first interior cell,
    corresponding to a zero-gradient (Neumann) upwind boundary condition.

    Parameters
    ----------
    array : ndarray
        Input array representing cell-centered values.

    Returns
    -------
    ndarray
        Array with one additional element prepended.
    """
    return np.insert(array, 0, array[0])


def churchill(eps: float, Dh: float, Re: np.ndarray) -> np.ndarray:
    """
    Compute the Darcy friction factor using the Churchill correlation.

    This correlation is valid for both laminar and turbulent flow
    and smoothly transitions between regimes.

    Parameters
    ----------
    eps : float
        Surface roughness.
    Dh : float
        Hydraulic diameter.
    Re : ndarray
        Reynolds number.

    Returns
    -------
    ndarray
        Darcy friction factor.
    """
    a = (2.457*np.log(1/((7/Re)**0.9 + (0.27*eps/Dh))))**16
    b = (37530/Re)**16

    f = 8*((8/Re)**12 + (1/((a+b)**1.5)))**(1/12)
    return f


def format_computation_time(duration: float) -> str:
    """
    Format a computation time in seconds into a human-readable string.

    Parameters
    ----------
    duration : float
        Elapsed CPU time in seconds.

    Returns
    -------
    str
        Formatted time string.
    """
    if duration >= 3600:
        hrs, r = divmod(duration, 3600)
        mins, secs = divmod(r, 60)
        return f"CPU Time: {hrs} [hr], {mins} [min], {int(secs)} [s]"

    elif duration >= 60:
        mins, secs = divmod(duration, 60)
        return f"CPU Time: {mins} [min], {int(secs)} [s]"

    else:
        return f"CPU Time: {duration:.1f} [s]"


def prandtl(k: np.ndarray, cp: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """
    Compute the Prandtl number.

    Parameters
    ----------
    k : float
        Thermal conductivity.
    cp : float
        Specific heat at constant pressure.
    mu : float
        Dynamic viscosity.

    Returns
    -------
    float
        Prandtl number.
    """
    return cp * mu / k


def residual(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> float:
    """
    Compute the normalized residual of a linear system.

    Parameters
    ----------
    A : numpy.ndarray
        System matrix.
    b : numpy.ndarray
        Right-hand-side vector.
    x : numpy.ndarray
        Solution vector.

    Returns
    -------
    float
        Relative residual norm.
    """
    return float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))


def reynolds(rho: np.ndarray, u: np.ndarray, L: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """
    Compute the Reynolds number.

    Parameters
    ----------
    rho : ndarray
        Fluid density.
    u : ndarray
        Flow velocity.
    L : ndarray
        Characteristic length.
    mu : ndarray
        Dynamic viscosity.

    Returns
    -------
    ndarray
        Reynolds number.
    """
    return rho * u * L / mu


def P_to_P0(P, u, rho):
    """
    Convert static pressure to stagnation pressure.

    Parameters
    ----------
    P : float
        Static pressure.
    u : float
        Flow velocity.
    rho : float
        Fluid density.

    Returns
    -------
    float
        Stagnation pressure.
    """
    return P + 0.5 * rho * u**2


def P0_to_P(P0: float, T: float, u: float, fluid: cedar.base.Fluid):
    """
    Convert stagnation pressure to static pressure via fixed-point iteration.

    Density is evaluated using the provided fluid model.

    Parameters
    ----------
    P0 : float
        Stagnation pressure.
    T : float
        Static temperature.
    u : float
        Flow velocity.
    fluid : cedar.base.Fluid
        Fluid property model.

    Returns
    -------
    float
        Static pressure.

    Notes
    -----
    Iteration terminates when convergence tolerance is met or
    an error is raised if convergence fails.
    """
    P = P0

    for _ in range(100):
        P_prev = P
        rho = fluid.rho_from_T_P(T, P)
        P = P0 - 0.5 * rho * u**2

        if abs(P - P_prev) <= 1:
            return P

    cedar.Log.error(
        f"P0_to_P did not converge: P0={P0:.0f}, T={T:.1f}, u={u:.1f}, "
        f"fluid={fluid.__class__.__name__}"
    )


def T_to_T0(T: float, u: float, cp: float):
    """
    Convert static temperature to stagnation temperature.

    Parameters
    ----------
    T : float
        Static temperature.
    u : float
        Flow velocity.
    cp : float
        Specific heat at constant pressure.

    Returns
    -------
    float
        Stagnation temperature.
    """
    return T + (0.5 * u**2) / cp


def T0_to_T(T0: float, P: float, u: float, fluid: cedar.base.Fluid):
    """
    Convert stagnation temperature to static temperature via iteration.

    Specific heat is evaluated using the provided fluid properties object.

    Parameters
    ----------
    T0 : float
        Stagnation temperature.
    P : float
        Static pressure.
    u : float
        Flow velocity.
    fluid : cedar.base.Fluid
        Fluid properties.

    Returns
    -------
    float
        Static temperature.
    """
    T = T0

    for _ in range(100):
        T_prev = T
        cp = fluid.cp_from_T_P(T, P)
        T = T0 - (0.5 * u**2) / cp

        if abs(T - T_prev) <= 1e-6:
            return T

    cedar.Log.error(
        f"T0_to_T did not converge: T0={T0:.0f}, P={P:.1f}, u={u:.1f}, "
        f"fluid={fluid.__class__.__name__}"
    )


def tetra_vol(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """
    Compute the volume of a tetrahedron.

    Parameters
    ----------
    p1, p2, p3, p4 : ndarray
        Vertex coordinates.

    Returns
    -------
    float
        Tetrahedron volume.
    """
    v1 = p1 - p4
    v2 = p2 - p4
    v3 = p3 - p4

    return np.abs(np.dot(v1, np.cross(v2, v3))) / 6.0


def tetra_center(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> np.ndarray:
    """
    Compute the centroid of a tetrahedron.

    Parameters
    ----------
    p1, p2, p3, p4 : ndarray
        Vertex coordinates.

    Returns
    -------
    ndarray
        Tetrahedron centroid.
    """
    return (p1 + p2 + p3 + p4) / 4


def triangle_area(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """
    Compute the area of a triangle.

    Parameters
    ----------
    p1, p2, p3 : ndarray
        Vertex coordinates.

    Returns
    -------
    float
        Triangle area.
    """
    v1 = p2 - p1
    v2 = p3 - p1

    return 0.5 * np.linalg.norm(np.cross(v1, v2))


def triangle_center(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """
    Compute the centroid of a triangle.

    Parameters
    ----------
    p1, p2, p3 : ndarray
        Vertex coordinates.

    Returns
    -------
    ndarray
        Triangle centroid.
    """
    return (p1 + p2 + p3) / 3


def triangle_normal(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """
    Compute the unit normal vector of a triangle.

    Parameters
    ----------
    p1, p2, p3 : ndarray
        Vertex coordinates.

    Returns
    -------
    ndarray
        Unit normal vector.
    """
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)


def westinghouse_modified_mccarthy_wolf(Re: np.ndarray, Pr: np.ndarray, T_wall: np.ndarray,
                                        T: np.ndarray, x: np.ndarray, d: np.ndarray) -> np.ndarray:
    """
    Compute the convective heat transfer coefficient using the
    modified McCarthy–Wolf correlation.

    Parameters
    ----------
    Re : ndarray
        Reynolds number.
    Pr : ndarray
        Prandtl number.
    T_wall : ndarray
        Wall temperature.
    T : ndarray
        Bulk fluid temperature.
    x : ndarray
        Axial distance.
    d : ndarray
        Characteristic diameter.

    Returns
    -------
    ndarray
        Nusselt-number-based heat transfer correlation value.
    """
    return 0.025 * (Re**0.8) * (Pr**0.4) * ((T_wall / T)**-0.55) * (1 + 0.3 * (x / d)**-0.7)