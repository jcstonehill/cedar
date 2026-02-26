import numpy as np
import meshio
import h5py as h5

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
        return f"{hrs} [hr], {mins} [min], {int(secs)} [s]"

    elif duration >= 60:
        mins, secs = divmod(duration, 60)
        return f"{mins} [min], {int(secs)} [s]"

    else:
        return f"{duration:.1f} [s]"

def gmsh_to_vtkhdf(gmsh_file: str, vtkhdf_file: str = None):
    VTK_TRIANGLE = 5
    VTK_TETRA = 10
    VERSION = (2,5)
    
    if vtkhdf_file is None:
        vtkhdf_file = gmsh_file.replace(".msh", ".vtkhdf")

    # Read the GMSH file using meshio
    mesh0: meshio.Mesh = meshio.read(gmsh_file)

    with h5.File(vtkhdf_file, "w") as f:
        root = f.create_group("VTKHDF", track_order = True)
        root.attrs["Version"] = VERSION
        root.attrs["Type"] = "MultiBlockDataSet"

        assembly = root.create_group("Assembly", track_order = True)

        # Delete useless GMSH specific data.
        initial_keys = list(mesh0.cell_sets.keys())
        for key in initial_keys:
            if key.startswith("gmsh:"):
                del mesh0.cell_sets[key]

        root.create_dataset("NumberOfPoints", data = (mesh0.points.shape[0],), dtype="i8")
        root.create_dataset("Points", data = mesh0.points, dtype="f")

        tri_pts = mesh0.cells_dict["triangle"]
        tet_pts = mesh0.cells_dict["tetra"]
        
        block_index = 0
        for cell_set_name in mesh0.cell_sets_dict:
            if "tetra" in mesh0.cell_sets_dict[cell_set_name]:
                block = root.create_group(cell_set_name, track_order = True)
                block.attrs["Version"] = VERSION
                block.attrs["Type"] = "UnstructuredGrid"
                block.attrs["Index"] = block_index
                block_index += 1

                block["NumberOfPoints"] = h5.SoftLink("/VTKHDF/NumberOfPoints")
                block["Points"] = h5.SoftLink("/VTKHDF/Points")
                
                mask = mesh0.cell_sets_dict[cell_set_name]["tetra"]
                block.create_dataset("NumberOfCells", data=(mask.shape[0],), dtype="i8")
                block.create_dataset("Types", data=np.full(mask.shape[0], VTK_TETRA), dtype="uint8")
                
                block.create_dataset("NumberOfConnectivityIds", data=(4*mask.shape[0],), dtype="i8")
                block.create_dataset("Connectivity", data=np.ravel(tet_pts[mask]), dtype="i8")
                block.create_dataset("Offsets", data = np.arange((mask.shape[0]+1)*4, step = 4), dtype="i8")

                assembly[cell_set_name] = h5.SoftLink("/VTKHDF/" + cell_set_name)

        for cell_set_name in mesh0.cell_sets_dict:
            if "triangle" in mesh0.cell_sets_dict[cell_set_name]:
                block = root.create_group(cell_set_name, track_order = True)
                block.attrs["Version"] = VERSION
                block.attrs["Type"] = "UnstructuredGrid"
                block.attrs["Index"] = block_index
                block_index += 1

                block["NumberOfPoints"] = h5.SoftLink("/VTKHDF/NumberOfPoints")
                block["Points"] = h5.SoftLink("/VTKHDF/Points")
                
                mask = mesh0.cell_sets_dict[cell_set_name]["triangle"]

                block.create_dataset("NumberOfCells", data=(mask.shape[0],), dtype="i8")
                block.create_dataset("Types", data=np.full(mask.shape[0], VTK_TRIANGLE), dtype="uint8")
                
                block.create_dataset("NumberOfConnectivityIds", data=(3*mask.shape[0],), dtype="i8")
                block.create_dataset("Connectivity", data=np.ravel(tri_pts[mask][:]), dtype="i8")
                block.create_dataset("Offsets", data = np.arange((mask.shape[0]+1)*3, step = 3), dtype="i8")

                assembly[cell_set_name] = h5.SoftLink("/VTKHDF/" + cell_set_name)

def MAPE(val, ref) -> float:
    """
    Compute the mean absolute percentage error (MAPE).

    Parameters
    ----------
    val : ndarray
        Computed or predicted values.
    ref : ndarray
        Reference values.

    Returns
    -------
    float
        Mean absolute percentage error in percent.

    Notes
    -----
    Zero values are replaced internally to avoid division by zero.
    """
    val = np.array(val)
    ref = np.array(ref)

    # Prevent division by zero
    ref[ref == 0] = 1e-64
    val[val == 0] = 1e-64

    return float(np.sum(np.abs((val - ref) / ref) * 100) / val.size)

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

def print_vtkhdf(file: str):
    try:
        with h5.File(file, 'r') as f:
            print(f"--- Contents of {file} ---")
            # Recursively print groups and datasets
            def print_h5_structure(name, obj):
                
                print(f"{'     ' * (obj.name.count('/')-1)}{name} ({type(obj).__name__})")

                if isinstance(obj, h5.Group):
                    for attr in obj.attrs:
                        print(f"{'     ' * (obj.name.count('/'))}{attr} : {obj.attrs[attr]}")

                    for key in obj:
                        print_h5_structure(key, obj[key])

                elif isinstance(obj, h5.Dataset):
                    print(f"{'     ' * (obj.name.count('/'))}Shape: {obj.shape}, Dtype: {obj.dtype}")

            # Start recursion from the root group
            for key in f:
                print_h5_structure(key, f[key])

    except FileNotFoundError:
        print(f"Error: File '{file}' not found.")
        
    except Exception as e:
        print(f"An error occurred: {e}")

def residual(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> float:
    """
    Compute the normalized residual of a linear system.

    Parameters
    ----------
    A : np.ndarray
        System matrix.
    b : np.ndarray
        Right-hand-side vector.
    x : np.ndarray
        Solution vector.

    Returns
    -------
    float
        Relative residual norm.
    """
    return np.linalg.norm(A @ x - b) / np.linalg.norm(b)


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


def P_to_P0(P: float, u: float, rho: float) -> float:
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


def P0_to_P(P0: float, u: float, rho: float) -> float:
    """
    Convert stagnation pressure to static pressure.

    Parameters
    ----------
    P0 : float
        Stagnation pressure.
    u : float
        Flow velocity.
    rho : float
        Fluid density.

    Returns
    -------
    float
        Static pressure.
    """

    return P0 - 0.5 * rho * u**2

def sort3(values):
    """
    Sort 3 numbers in ascending order (computationally efficient).

    """
    if len(values) != 3:
        raise Exception("sort3 can only handle an iterable of size 3.")
    
    a = values[0]
    b = values[1]
    c = values[2]

    if a > b: a, b = b, a
    if b > c: b, c = c, b
    if a > b: a, b = b, a
    
    return (a, b, c)

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


def T0_to_T(T0: float, u: float, cp: float):
    """
    Convert stagnation temperature to static.

    Parameters
    ----------
    T0 : float
        Stagnation temperature.
    u : float
        Flow velocity.
    cp : float
        Specific heat capacity.

    Returns
    -------
    float
        Static temperature.
    """
    return T0 - (0.5 * u**2) / cp


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


# def tetra_center(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> np.ndarray:
#     """
#     Compute the centroid of a tetrahedron.

#     Parameters
#     ----------
#     p1, p2, p3, p4 : ndarray
#         Vertex coordinates.

#     Returns
#     -------
#     ndarray
#         Tetrahedron centroid.
#     """
#     return (p1 + p2 + p3 + p4) / 4


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


# def westinghouse_modified_mccarthy_wolf(Re: np.ndarray, Pr: np.ndarray, T_wall: np.ndarray,
#                                         T: np.ndarray, x: np.ndarray, d: np.ndarray) -> np.ndarray:
#     """
#     Compute the convective heat transfer coefficient using the
#     modified McCarthy–Wolf correlation.

#     Parameters
#     ----------
#     Re : ndarray
#         Reynolds number.
#     Pr : ndarray
#         Prandtl number.
#     T_wall : ndarray
#         Wall temperature.
#     T : ndarray
#         Bulk fluid temperature.
#     x : ndarray
#         Axial distance.
#     d : ndarray
#         Characteristic diameter.

#     Returns
#     -------
#     ndarray
#         Nusselt-number-based heat transfer correlation value.
#     """
#     return 0.025 * (Re**0.8) * (Pr**0.4) * ((T_wall / T)**-0.55) * (1 + 0.3 * (x / d)**-0.7)