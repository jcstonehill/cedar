from abc import ABC, abstractmethod
import h5py as h5
import itertools
import numpy as np
import os
import sys
import itertools

import cedar


# MESH SHOULD USE A CUSTOM ID TO DIFFERENTIATE POINT SETS IN THE VTKHDF FILES RATHER THAN THE MODEL NAME.
# CAN I IMPLEMENT THIS TO BE AUTOMATIC?

class Mesh(ABC):

    id: itertools.count = itertools.count()

    def __init__(self):
        self.id: int = next(Mesh.id)
        self.dim: int = 0

        self.regions: list[str] = []
        self.boundaries: list[str] = []

        self.N_pts: np.ndarray = 0
        self.pts: np.ndarray = []

        self.N_cells: int = 0
        self.cell_pts_i: np.ndarray = []
        self.cell_faces_i: np.ndarray = []

        self.N_faces: int = 0
        self.face_pts_i: np.ndarray = []
        self.face_cells_i: np.ndarray = []
        self.face_is_interior: np.ndarray = []

        self.face_n: np.ndarray = []
        self.face_d: np.ndarray = []
        self.face_w: np.ndarray = []

        self.cell_centers: np.ndarray = []
        self.face_centers: np.ndarray = []
        
        self.region_i: dict[str, np.ndarray] = {}
        self.region_N: dict[str, int] = {}

        self.boundary_i: dict[str, np.ndarray] = {}
        self.boundary_N: dict[str, int] = {}

    def get_local_cell_i(self, global_cell_i: int) -> tuple[str, int]:
        for region in self.regions:
            cell_indices = self.region_i[region]

            if global_cell_i in cell_indices:
                return region, np.where(cell_indicies = global_cell_i)[0][0]
        
        raise Exception("Can't find this global_cell_i in any region: " + str(global_cell_i))

    @abstractmethod
    def build(self):
        pass

    @abstractmethod
    def vtkhdf_dict(self) -> dict:
        pass

class Mesh0D(Mesh):
    def __init__(self):
        self.id: int = next(Mesh.id)
        self.dim: int = 0

        raise NotImplementedError()

class Mesh1D(Mesh):
    basis_vec = {
        "x"  : np.array([ 1,  0,  0]),
        "y"  : np.array([ 0,  1,  0]),
        "z"  : np.array([ 0,  0,  1]),
        "-x" : np.array([-1,  0,  0]),
        "-y" : np.array([ 0, -1,  0]),
        "-z" : np.array([ 0,  0, -1])
    }

    def __init__(
            self,
            N_cells: int,
            L: float,
            region: str,
            start_boundary: str,
            end_boundary: str,
            start: tuple[float],
            basis: str = "z"
        ):
        self.id: int = next(Mesh.id)
        self.dim: int = 1

        self.L = L
        self.start = np.array(start)
        self.end = start + L*self.basis_vec[basis]
        self.basis = basis

        self.ds = L/N_cells

        self.regions: list[str] = [region]
        self.boundaries: list[str] = [start_boundary, end_boundary]

        # Attributes After Build
        self.N_pts: np.ndarray = 0
        self.pts: np.ndarray = []

        self.N_cells: int = N_cells
        self.cell_pts_i: np.ndarray = []
        self.cell_faces_i: np.ndarray = []

        self.N_faces: int = 0
        self.face_pts_i: np.ndarray = []
        self.face_cells_i: np.ndarray = []
        self.face_is_interior: np.ndarray = []

        self.face_n: np.ndarray = []
        self.face_d: np.ndarray = []
        self.face_w: np.ndarray = []

        self.cell_centers: np.ndarray = []
        self.face_centers: np.ndarray = []
        
        self.region_i: dict[str, np.ndarray] = {}
        self.region_N: dict[str, int] = {}

        self.boundary_i: dict[str, np.ndarray] = {}
        self.boundary_N: dict[str, int] = {}

    def build(self):
        # Attributes After Build
        self.N_pts: np.ndarray = self.N_cells + 1

        delta = self.end - self.start
        self.pts = np.zeros((self.N_pts, 3), dtype = np.float64)
        for i in range(self.N_pts):
            self.pts[i][:] = self.start + delta*i/(self.N_pts-1)

        self.N_cells: int = self.N_cells # Already provided by the user.
        self.cell_pts_i: np.ndarray = np.array([(i, i+1) for i in range(self.N_cells)])
        self.cell_faces_i: np.ndarray = np.copy(self.cell_pts_i) # For 1D, points and faces are the same thing

        # For 1D, points and faces are the same thing
        self.N_faces: int = self.N_pts          
        self.face_pts_i: np.ndarray = np.arange(self.N_pts)
        self.face_is_interior: np.ndarray = np.full(self.N_faces, True)
        self.face_is_interior[0] = False
        self.face_is_interior[-1] = False
        self.face_cells_i = np.zeros((self.N_faces, 2), dtype=np.int64)
        self.face_cells_i[-1,:] = [self.N_cells-1, self.N_cells-1]

        self.cell_centers: np.ndarray = np.zeros((self.N_cells, 3), dtype = np.float64)
        delta = delta / np.linalg.norm(delta)
        for i in range(self.N_cells):
            self.cell_centers[i][:] = self.start + delta*self.ds*(0.5+i)
        self.face_centers: np.ndarray = np.copy(self.pts) # For 1D, points and faces are the same thing

        self.face_n: np.ndarray = np.full((self.N_faces, 3), self.basis_vec[self.basis], dtype = np.float64)
        self.face_d: np.ndarray = np.full((self.N_faces, 2), self.ds/2, dtype = np.float64)
        self.face_d[0][1] = 0
        self.face_d[-1][1] = 0
        self.face_w: np.ndarray = np.full(self.N_faces, 0.5)
        self.face_w[0] = 1
        self.face_w[-1] = 1

        self.region_i: dict[str, np.ndarray] = {self.regions[0] : np.arange(self.N_cells)}
        self.region_N: dict[str, int] = {self.regions[0] : self.N_cells}

        self.boundary_i: dict[str, np.ndarray] = {
            self.boundaries[0] : np.array([0]),
            self.boundaries[1] : np.array([-1])
        }
        self.boundary_N: dict[str, int] = {
            self.boundaries[0] : 1,
            self.boundaries[1] : 1
        }

    def vtkhdf_dict(self) -> dict:
        vtkhdf_dict = {
            f"{self.id}_NumberOfPoints" : (self.pts.shape[0],),
            f"{self.id}_Points" : self.pts,
            "blocks" : {}
        }

        vtkhdf_dict["blocks"][self.regions[0]] = {
            "NumberOfCells" : ((self.N_cells,), "i8"),
            "Types" : (np.full(self.N_cells, 3), "uint8"), # 3 means VTK_LINE
            "NumberOfConnectivityIds" : ((2*self.N_cells,), "i8"),
            "Connectivity" : (np.ravel(self.cell_pts_i), "i8"),
            "Offsets" : (np.arange(self.N_pts*2, step = 2), "i8")
        }

        vtkhdf_dict["blocks"][self.boundaries[0]] = {
            "NumberOfCells" : ((1,), "i8"),
            "Types" : ((1,), "uint8"), # 1 means VTK_VERTEX
            "NumberOfConnectivityIds" : ((1,), "i8"),
            "Connectivity" : ((0,), "i8"),
            "Offsets" : ((0,1), "i8")
        }

        vtkhdf_dict["blocks"][self.boundaries[1]] = {
            "NumberOfCells" : ((1,), "i8"),
            "Types" : ((1,), "uint8"), # 1 means VTK_VERTEX
            "NumberOfConnectivityIds" : ((1,), "i8"),
            "Connectivity" : ((self.N_cells,), "i8"),
            "Offsets" : ((0,1), "i8")
        }

        return vtkhdf_dict


class Mesh3D(Mesh):
    def __init__(self, file: str):
        self.id: int = next(Mesh.id)
        self.dim: int = 3

        self.file: str = file

        self.regions: list[str] = self.get_regions()
        self.boundaries: list[str] = self.get_boundaries()

        # Attributes After Build
        self.N_pts: np.ndarray = 0
        self.pts: np.ndarray = []

        self.N_cells: int = 0
        self.cell_pts_i: np.ndarray = []
        self.cell_faces_i: np.ndarray = []

        self.N_faces: int = 0
        self.face_pts_i: np.ndarray = []
        self.face_cells_i: np.ndarray = []
        self.face_is_interior: np.ndarray = []

        self.face_n: np.ndarray = []
        self.face_d: np.ndarray = []
        self.face_w: np.ndarray = []

        self.cell_centers: np.ndarray = []
        self.face_centers: np.ndarray = []

        self.region_i: dict[str, np.ndarray] = {}
        self.region_N: dict[str, int] = {}

        self.boundary_i: dict[str, np.ndarray] = {}
        self.boundary_N: dict[str, int] = {}

        # Specific to 3D
        self.cell_vols: np.ndarray = []
        self.face_areas: np.ndarray = []
        self.region_vol: dict[str, float] = {}
        self.boundary_area: dict[str, float] = {}

    def build(self):
        with h5.File(self.file, "r") as f:
            root = f["VTKHDF"]

            self.N_pts = root["NumberOfPoints"][:]
            self.pts = root["Points"][:]
            
            face_dict = self._create_cells(root)
            self._create_faces(face_dict)
            self._define_boundaries(root, face_dict)
            self._calc_cell_data()
            self._calc_face_data()

    def get_boundaries(self) -> list[str]:
        boundaries: list[str] = []

        with h5.File(self.file, "r") as f:
            root = f["VTKHDF"]

            for block, data in root["Assembly"].items():
                cell_type = data["Types"][0]

                if not np.all(cell_type == data["Types"][:]):
                    raise Exception("A block can only have one cell type.")

                elif cell_type == 5: # VTK_TRIANGLE
                    boundaries.append(block)
        
        return boundaries
    
    def get_regions(self) -> list[str]:
        regions: list[str] = []

        with h5.File(self.file, "r") as f:
            root = f["VTKHDF"]

            for block, data in root["Assembly"].items():
                cell_type = data["Types"][0]

                if not np.all(cell_type == data["Types"][:]):
                    raise Exception("A block can only have one cell type.")

                if cell_type == 10: # VTK_TETRA
                    regions.append(block)
                
        return regions

    def _calc_cell_data(self):
        # 1. Fetch all vertex coordinates at once using 'fancy indexing'
        # Shape becomes: (N_cells, 4 vertices, 3 coords)
        cells = self.pts[self.cell_pts_i]

        # 2. Calculate vectors from the 4th vertex (p4) to the others (p1, p2, p3)
        # We subtract the 4th vertex (index 3) from the first three (indices 0, 1, 2)
        v1 = cells[:, 0, :] - cells[:, 3, :]
        v2 = cells[:, 1, :] - cells[:, 3, :]
        v3 = cells[:, 2, :] - cells[:, 3, :]

        # 3. Compute Cross Product (v2 x v3) for all cells
        cross_prod = np.cross(v2, v3, axis=1)

        # 4. Compute Dot Product (v1 . (v2 x v3)) for all cells
        # We use np.einsum for a fast batch dot product, or simple element-wise mult + sum
        scalar_triple = np.sum(v1 * cross_prod, axis=1)

        # 5. Final Volume Calculation
        self.cell_vols = np.abs(scalar_triple) / 6.0
        
        # 6. Vectorized calculation for cell centers (mean of the 4 vertices)
        self.cell_centers = np.mean(cells, axis=1)

        for region in self.regions:
            self.region_vol[region] = np.sum(self.cell_vols[self.region_i[region]])

    def _calc_face_data(self):
        # 1. Gather all points at once using "Fancy Indexing"
        # p1, p2, p3 will be arrays of shape (N_faces, 3)
        p1 = self.pts[self.face_pts_i[:, 0]]
        p2 = self.pts[self.face_pts_i[:, 1]]
        p3 = self.pts[self.face_pts_i[:, 2]]

        # 2. Vectorized Centers
        # Compute the mean of the coordinates for all faces simultaneously
        self.face_centers = (p1 + p2 + p3) / 3.0

        # 3. Vectorized Normals and Areas
        # Calculate edge vectors for every face
        edge1 = p2 - p1
        edge2 = p3 - p1

        # Cross product of all faces at once (N, 3)
        cross_prod = np.cross(edge1, edge2)

        # The magnitude (norm) of the cross product is 2 * Area
        # axis=1 calculates the norm across the [x, y, z] vector for each row
        cross_norms = np.linalg.norm(cross_prod, axis=1)

        self.face_areas = 0.5 * cross_norms.astype(np.float64)

        # Normalize the vectors to get unit normals
        # We add [:, np.newaxis] to allow division of (N,3) by (N,1)
        # Added epsilon (1e-12) to prevent division by zero for degenerate faces
        self.face_n = (cross_prod / (cross_norms[:, np.newaxis] + 1e-12)).astype(np.float64)

        # 4. Vectorized Neighbor Distances (face_d)
        # Gather cell centers for all adjacent cells
        c1_centers = self.cell_centers[self.face_cells_i[:, 0]]
        c2_centers = self.cell_centers[self.face_cells_i[:, 1]]

        # Calculate vectors from face centers to cell centers
        v1 = c1_centers - self.face_centers
        v2 = c2_centers - self.face_centers

        # Since np.dot is for matrix multiplication, we use element-wise multiply + sum 
        # to get dot products for rows of vectors.
        d1 = np.abs(np.sum(v1 * self.face_n, axis=1))
        d2 = np.abs(np.sum(v2 * self.face_n, axis=1))
        d2[~self.face_is_interior] = 0

        self.face_d = np.stack((d1, d2), axis=1).astype(np.float64)

        self.face_w = 1 / (1+d1/d2)

        for boundary in self.boundaries:
            self.boundary_area[boundary] = np.sum(self.face_areas[self.boundary_i[boundary]])

    def _create_cells(self, h5_root: h5.Group) -> dict:
        for region in self.regions:
            data = h5_root["Assembly"][region]

            N = data["NumberOfCells"][0]
            cell_pts_i = np.reshape(data["Connectivity"][:], (N, data["Offsets"][1]))

            self.region_i[region] = np.arange(self.N_cells, self.N_cells+N)
            self.N_cells += N
            self.region_N[region] = N

            self.cell_pts_i.extend(cell_pts_i)
        
        # Now we know N_cells, so we can allocate data arrays.
        self.cell_pts_i: np.ndarray = np.array(self.cell_pts_i, dtype=np.uint32)
        self.cell_faces_i: np.ndarray = np.zeros((self.N_cells, 4), dtype = np.uint32)
        self.cell_centers: np.ndarray = np.zeros((self.N_cells, 3), dtype = np.float64)

        # Finally, we need to create a dictionary of unique faces. key is tuple
        # of (p1, p2, p3) where p1, p2, p3 are face point i in ascending order.
        # face_dict[key][0] = face_i
        # face_dict[key][1] = face_pts_i
        # face_dict[key][2] = cell_1_i
        # face_dict[key][3] = cell_2_i

        face_dict = {}
        new_face_i = 0

        for cell_i, cell_pt_set in enumerate(self.cell_pts_i):
            combinations = itertools.combinations(cell_pt_set, 3)

            for i, face_pt_set in enumerate(combinations):
                face_key = cedar.functions.sort3(face_pt_set)

                if face_key not in face_dict:
                    face_i = new_face_i
                    new_face_i += 1

                    face_dict[face_key] = [face_i, face_pt_set, cell_i, None]

                else:
                    face_i = face_dict[face_key][0]
                    face_dict[face_key][3] = cell_i

                self.cell_faces_i[cell_i][i]

        return face_dict

    def _create_faces(self, face_dict: dict):
        self.N_faces = len(face_dict)
        self.face_pts_i = np.zeros((self.N_faces, 3), dtype = np.uint32)
        self.face_cells_i = np.zeros((self.N_faces, 2), dtype = np.uint32)
        self.face_is_interior = np.full(self.N_faces, False, dtype = np.bool)

        for face_i, face_pts_i, cell1_i, cell2_i in face_dict.values():
            self.face_pts_i[face_i] = face_pts_i

            self.face_cells_i[face_i][0] = cell1_i
            self.face_cells_i[face_i][1] = cell2_i if cell2_i is not None else cell1_i

            self.face_is_interior[face_i] = True if cell2_i is not None else False

    def _define_boundaries(self, h5_root: h5.Group, face_dict: dict):
        for boundary in self.boundaries:
            # TODO Can I just use the boundary key in the root to get all of the
            # data instead of going through Assembly?
            data = h5_root["Assembly"][boundary]

            N = data["NumberOfCells"][0]
            face_pts_i = np.reshape(data["Connectivity"][:], (N, data["Offsets"][1]))

            self.boundary_i[boundary] = np.zeros(N, dtype = np.uint32)
            self.boundary_N[boundary] = N

            for i, pt_set in enumerate(face_pts_i):
                key = cedar.functions.sort3(pt_set)
                self.boundary_i[boundary][i] = face_dict[key][0]

    def vtkhdf_dict(self):
        vtkhdf_dict = {
            f"{self.id}_NumberOfPoints" : (self.pts.shape[0],),
            f"{self.id}_Points" : self.pts,
            "blocks" : {}
        }

        for region in self.regions:
            vtkhdf_dict["blocks"][region] = {
                "NumberOfCells" : ((self.region_N[region],), "i8"),
                "Types" : (np.full(self.region_N[region], 10), "uint8"), # 10 means VTK_TETRA
                "NumberOfConnectivityIds" : ((4*self.region_N[region],), "i8"),
                "Connectivity" : (np.ravel(self.cell_pts_i[self.region_i[region]]), "i8"),
                "Offsets" : (np.arange((self.region_N[region]+1)*4, step = 4), "i8")
            }

        for boundary in self.boundaries:
            vtkhdf_dict["blocks"][boundary] = {
                "NumberOfCells" : ((self.boundary_N[boundary],), "i8"),
                "Types" : (np.full(self.boundary_N[boundary], 5), "uint8"), # 5 means VTK_TRIANGLE
                "NumberOfConnectivityIds" : ((3*self.boundary_N[boundary],), "i8"),
                "Connectivity" : (np.ravel(self.face_pts_i[self.boundary_i[boundary]]), "i8"),
                "Offsets" : (np.arange((self.boundary_N[boundary]+1)*3, step = 3), "i8")
            }

        return vtkhdf_dict