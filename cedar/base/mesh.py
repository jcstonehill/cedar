from abc import ABC, abstractmethod
import numpy as np


class Mesh(ABC):
    """
    Abstract base class for Meshes.

    This class serves as a common interface for concrete mesh
    implementations (e.g., 1D, 2D, or 3D meshes). It stores mesh geometry,
    cell/face connectivity, and optional region and boundary groupings.
    No assumptions are made about dimensionality beyond what is encoded
    in the stored arrays.

    Concrete subclasses are responsible for populating all attributes
    consistently.

    Attributes
    ----------
    points : numpy.ndarray
        Array of mesh point coordinates with shape ``(N_points, 3)``.

    N_cells : int
        Number of cells in the mesh.

    cell_point_ids : numpy.ndarray
        Integer array mapping each cell to its defining point indices.

    cell_face_ids : numpy.ndarray
        Integer array mapping each cell to its associated face indices.

    cell_centers : numpy.ndarray
        Coordinates of cell centers with shape ``(N_cells, 3)``.

    N_faces : int
        Number of faces in the mesh.

    face_point_ids : numpy.ndarray
        Integer array mapping each face to its defining point indices.

    face_cell_ids : numpy.ndarray
        Integer array mapping each face to the adjacent cell indices.

    face_is_interior : numpy.ndarray
        Boolean array indicating whether each face is an interior face.

    face_centers : numpy.ndarray
        Coordinates of face centers.

    face_n : numpy.ndarray
        Face normal vectors.

    face_L : numpy.ndarray
        Characteristic face length(s).

    face_w : numpy.ndarray
        Face weighting or interpolation coefficients.

    regions : dict
        Mapping of region names to region metadata.

    region_cell_ids : dict
        Mapping of region names to arrays of cell indices.

    region_N : dict
        Mapping of region names to the number of cells in each region.

    boundaries : dict
        Mapping of boundary names to boundary metadata.

    boundary_face_ids : dict
        Mapping of boundary names to arrays of face indices.

    boundary_N : dict
        Mapping of boundary names to the number of faces on each boundary.
    """

    def __init__(self):
        self.points: np.ndarray = None

        self.N_cells: np.ndarray = None
        self.cell_point_ids: np.ndarray = None
        self.cell_face_ids: np.ndarray = None
        
        self.cell_centers: np.ndarray = None

        self.N_faces: np.ndarray = None
        self.face_point_ids: np.ndarray = None
        self.face_cell_ids: np.ndarray = None
        self.face_is_interior: np.ndarray = None

        self.face_centers: np.ndarray = None
        self.face_n = None
        self.face_L = None
        self.face_w = None

        self.regions = {}
        self.region_cell_ids = {}
        self.region_N = {}

        self.boundaries = {}
        self.boundary_face_ids = {}
        self.boundary_N = {}