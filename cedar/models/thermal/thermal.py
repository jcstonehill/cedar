import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import cg, spilu, LinearOperator, aslinearoperator, qmr, bicgstab
from scipy.linalg import solve as scipy_solve
from dataclasses import dataclass

import cedar


@dataclass
class Params:
    materials: dict[cedar.base.Material]
    bc_types: dict[int]

@dataclass
class Vars:
    Qdots: dict[cedar.ScalarVar]
    Qdot_shapes: dict[cedar.Mesh3DVar]
    bcs: dict[cedar.Mesh3DVar]

    T: cedar.Mesh3DVar
    vol_Qdot: cedar.Mesh3DVar
    k: cedar.Mesh3DVar
    rho: cedar.Mesh3DVar
    cp: cedar.Mesh3DVar

    boundary_Qdots: dict[cedar.Mesh3DVar]
    boundary_Ts: dict[cedar.Mesh3DVar]
    
class Thermal(cedar.base.Model):
    """
    Solves the heat conduction equation.

    Used to calculate temperature distributions in solids.
    """

    model_type = "Thermal"

    bc_info_dict = {
        "Dirichlet" : (1, "K"),
        "dirichlet" : (1, "K"),
        "Fixed Value" : (1, "K"),
        "fixed value" : (1, "K"),

        "Neumann" : (2, "W/m2"),
        "neumann" : (2, "W/m2"),
        "Flux" : (2, "W/m2"),
        "flux" : (2, "W/m2"),
    }

    @property
    def mesh(self) -> cedar.Mesh3D:
        return self._mesh
    
    @property
    def params(self) -> Params:
        return self._params
    
    @property
    def vars(self) -> Vars:
        return self._vars
    
    def __init__(self, name, mesh: cedar.Mesh3D, materials: dict[cedar.base.Material] = {}):
        self.name = name
        self._mesh = mesh

        Qdots = {}
        Qdot_shapes = {}
        bc_types = {}
        bcs = {}

        boundary_Qdots = {}
        boundary_Ts = {}

        for region in self.mesh.regions:
            if region not in materials.keys():
                materials[region] = None

            Qdots[region] = self.add_var(cedar.ScalarVar(self.name, "Qdot_" + region, 0, "W", self.mesh, region))
            Qdot_shapes[region] = self.add_var(cedar.Mesh3DVar(self.name, "Qdot_shape_" + region, 1, mesh = self.mesh, region = region))

        for boundary in self.mesh.boundaries:
            bc_types[boundary] = "0"
            bcs[boundary] = self.add_var(cedar.Mesh3DVar(self.name, "bc_" + boundary, 0, "W/m2", mesh = self.mesh, boundary = boundary))
            boundary_Qdots[boundary] = self.add_var(cedar.Mesh3DVar(self.name, "boundary_Qdot_" + boundary, 0, "W", self.mesh, boundary = boundary))
            boundary_Ts[boundary] = self.add_var(cedar.Mesh3DVar(self.name, "boundary_T_" + boundary, 300, "K", self.mesh, boundary = boundary))

        self._params = Params(
            materials = materials,
            bc_types = bc_types
        )

        self._vars = Vars(
            Qdots = Qdots,
            Qdot_shapes = Qdot_shapes,
            bcs = bcs,

            T = self.add_var(cedar.Mesh3DVar(self.name, "T", 300, "K", self.mesh)),
            vol_Qdot = self.add_var(cedar.Mesh3DVar(self.name, "vol_Qdot", 0, "W/m3", self.mesh)),
            k = self.add_var(cedar.Mesh3DVar(self.name, "k", 1, "W/m-K", self.mesh)),
            rho = self.add_var(cedar.Mesh3DVar(self.name, "rho", 1, "kg/m3", self.mesh)),
            cp = self.add_var(cedar.Mesh3DVar(self.name, "cp", 1, "J/kg-K", self.mesh)),
            boundary_Qdots = boundary_Qdots,
            boundary_Ts = boundary_Ts
        )

    def check(self):
        for region in self.mesh.regions.keys():
            if region not in self.params.materials.keys():
                self.log_error("No material assigned to region: \"" + str(region) + "\"")

    def cp_by_cell(self, T):
        cp = np.empty(self.mesh.N_cells)

        for region, mask in self.mesh.regions.items():
            mat = self.params.materials[region]
            cp[mask] = mat.cp(T[mask])

        return cp

    def k_by_cell(self, T):
        k = np.empty(self.mesh.N_cells)
        
        for region, mask in self.mesh.regions.items():
            mat = self.params.materials[region]
            k[mask] = mat.k(T[mask])

        return k

    def k_by_face(self, T):
        cell_k = self.k_by_cell(T)
        k = np.empty(self.mesh.N_faces)

        f = self.mesh.face_is_interior
        c1 = self.mesh.face_cell_ids[f, 0]
        c2 = self.mesh.face_cell_ids[f, 1]
        w  = self.mesh.face_w[f]

        k[f] = w * cell_k[c1] + (1 - w) * cell_k[c2]

        # boundary faces
        bf = ~f
        k[bf] = cell_k[self.mesh.face_cell_ids[bf, 0]]

        return k

    def rho_by_cell(self):
        rho = np.empty(self.mesh.N_cells)
        
        for region, mask in self.mesh.regions.items():
            mat = self.params.materials[region]
            rho[mask] = mat.rho_rt()

        return rho
    
    def S(self):
        Qdot = np.zeros(self.mesh.N_cells)

        for region in self.vars.Qdots:
            shape = self.vars.Qdot_shapes[region].val
            w = self.mesh.cell_vols[self.mesh.regions[region]]*shape
            w = w / sum(w)

            Qdot[self.mesh.regions[region]] = self.vars.Qdots[region].val*w

        return Qdot

    def set_bc(self, boundary, bc_type, bc_val):
        bc_type, units = self.bc_info_dict[bc_type]

        self.params.bc_types[boundary] = bc_type
        self.vars.bcs[boundary].set(bc_val)
        self.vars.bcs[boundary].units = units

    def set_Qdot(self, Qdot, region = None):
            if region is None:
                for region_name in self.mesh.regions:
                    self.vars.Qdots[region_name].set(Qdot*self.mesh.region_vols[region_name]/self.mesh.vol)

            else:
                self.vars.Qdots[region_name].set(Qdot)
    
    def set_material(self, region, material):
        self.params.materials[region] = material

    def setup(self):
        self.vars.T.set(self.vars.T.initial)

    def iterate(self, t0 = None, dt = None, res_reduc = 1e-2):
        # relax = 0.2

        T_initial = self.vars.T.initial
        T = self.vars.T.val

        if(dt is not None):
            rho = self.rho_by_cell()
            mass = rho*self.mesh.cell_vols

        S = self.S()

        face_k = self.k_by_face(T)

        if(dt is None):
            A, b = self._steady_matrices_sparse(face_k, S)

        else:
            cp = self.cp_by_cell(T)
            A, b = self._unsteady_matrices_sparse(dt, mass, face_k, cp, T_initial, S)
        
        res = cedar.helper.residual(A, b, T)
        tgt_res = res_reduc * res
        #tgt_res = 1e-8
        T, _ = cg(A, b, T, rtol = tgt_res)

        if _ != 0:
            self.log_error("Thermal model did not converge.")

        self._post(T, S/self.mesh.cell_vols, self.rho_by_cell(), self.cp_by_cell(T), self.k_by_cell(T), self.k_by_face(T))
        return tgt_res
            # #print(res, tgt_res)
            # if res <= tgt_res:
            #     #cedar.Log.message("Converged!")
                
            #     return res

            

            # face_k = self.k_by_face(T)

            # if(dt is None):
            #     A, b = self._steady_matrices_sparse(face_k, S)

            # else:
            #     cp = self.cp_by_cell(T)
            #     A, b = self._unsteady_matrices_sparse(dt, mass, face_k, cp, T_initial, S)

            # res = cedar.helper.residual(A, b, T)

            #cedar.Log.message(f"{i+1:>3} |R|: {res:.4e}")
            
            
            #T, _ = cg(A, b, T, rtol = 1e-12)
            #print(T, _)

    def _post(self, T, vol_Qdot, rho, cp, k, face_k):
        for boundary in self.mesh.boundaries:

            face_ids = self.mesh.boundaries[boundary]

            face_data = zip(
                self.mesh.face_cell_ids[face_ids],
                face_k[face_ids],
                self.mesh.face_areas[face_ids],
                self.mesh.face_L[face_ids]
            )

            bc_type = self.params.bc_types[boundary]
            bc_vals = self.vars.bcs[boundary].val

            boundary_T = np.zeros(self.mesh.boundary_N[boundary])
            boundary_Qdot = np.zeros(self.mesh.boundary_N[boundary])

            for i, (cell_ids, k_i, area, L) in enumerate(face_data):
                C = k_i * area / L
                c1 = cell_ids[0]

                # Adiabatic
                if bc_type == 0:
                    boundary_T[i] = T[c1]
                    boundary_Qdot[i] = 0

                # Dirichlet
                if bc_type == 1:
                    boundary_T[i] = bc_vals[i]
                    boundary_Qdot[i] = C*(T[c1] - bc_vals[i])
                    
                # Neumann
                if bc_type == 2:
                    boundary_T[i] = T[c1] - bc_vals[i]*area/C
                    boundary_Qdot[i] = bc_vals[i]*area

            self.vars.boundary_Qdots[boundary].set(boundary_Qdot)
            self.vars.boundary_Ts[boundary].set(boundary_T)

        self.vars.T.set(T)
        self.vars.vol_Qdot.set(vol_Qdot)
        self.vars.rho.set(rho)
        self.vars.cp.set(cp)
        self.vars.k.set(k)

    def _steady_matrices_sparse(self, face_k, S):
        N = self.mesh.N_cells
        b = np.copy(S)

        rows = []
        cols = []
        data = []

        # -----------------------------------------
        #  INTERIOR FACES
        # -----------------------------------------
        face_data = zip(
            self.mesh.face_cell_ids[self.mesh.face_is_interior],
            face_k[self.mesh.face_is_interior],
            self.mesh.face_areas[self.mesh.face_is_interior],
            self.mesh.face_L[self.mesh.face_is_interior]
        )

        for cell_ids, k, area, L in face_data:
            Df = k * area / L
            c1, c2 = cell_ids

            # A[c1][c1] += C
            rows.append(c1); cols.append(c1); data.append( Df)
            # A[c1][c2] += -C
            rows.append(c1); cols.append(c2); data.append(-Df)
            # A[c2][c2] += C
            rows.append(c2); cols.append(c2); data.append( Df)
            # A[c2][c1] += -C
            rows.append(c2); cols.append(c1); data.append(-Df)

        # -----------------------------------------
        #  BOUNDARY FACES
        # -----------------------------------------
        for boundary in self.mesh.boundaries:

            face_ids = self.mesh.boundaries[boundary]

            face_data = zip(
                self.mesh.face_cell_ids[face_ids],
                face_k[face_ids],
                self.mesh.face_areas[face_ids],
                self.mesh.face_L[face_ids]
            )

            bc_type = self.params.bc_types[boundary]
            bc_vals = self.vars.bcs[boundary].val

            for i, (cell_ids, k, area, L) in enumerate(face_data):
                c1 = cell_ids[0]

                # Dirichlet
                if bc_type == 1:  
                    Df = k * area / L

                    rows.append(c1); cols.append(c1); data.append(Df)
                    b[c1] += Df * bc_vals[i]

                # Neumann
                elif bc_type == 2:
                    b[c1] -= bc_vals[i] * area

        # -----------------------------------------
        # Build sparse matrix
        # -----------------------------------------
        A = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

        return A, b

    def _unsteady_matrices_sparse(self, dt, mass, face_k, cp, T_initial, S):
        # lists for COO sparse matrix storage
        rows = []
        cols = []
        data = []

        N = self.mesh.N_cells
        b = np.zeros(N)

        # ------------------------------
        #   Time + source terms
        # ------------------------------
        for i in range(N):
            diag_val = mass[i] * cp[i] / dt
            rows.append(i); cols.append(i); data.append(diag_val)

            b[i] += S[i] + T_initial[i] * diag_val

        # ------------------------------
        #   INTERIOR FACES
        # ------------------------------
        face_data = zip(self.mesh.face_cell_ids[self.mesh.face_is_interior],
                        face_k[self.mesh.face_is_interior],
                        self.mesh.face_areas[self.mesh.face_is_interior],
                        self.mesh.face_L[self.mesh.face_is_interior])

        for cell_ids, k, area, L in face_data:
            c1, c2 = cell_ids
            Df = k * area / L

            # A[c1][c1] += Df
            rows.append(c1); cols.append(c1); data.append(Df)
            # A[c2][c2] += Df
            rows.append(c2); cols.append(c2); data.append(Df)

            # A[c1][c2] -= Df
            rows.append(c1); cols.append(c2); data.append(-Df)
            # A[c2][c1] -= Df
            rows.append(c2); cols.append(c1); data.append(-Df)

        # ------------------------------
        #   BOUNDARY FACES
        # ------------------------------
        for boundary in self.mesh.boundaries:

            face_data = zip(self.mesh.face_cell_ids[self.mesh.boundaries[boundary]],
                            face_k[self.mesh.boundaries[boundary]],
                            self.mesh.face_areas[self.mesh.boundaries[boundary]],
                            self.mesh.face_L[self.mesh.boundaries[boundary]])

            for i, (cell_ids, k, area, L) in enumerate(face_data):
                c1 = cell_ids[0]

                # Dirichlet
                if self.params.bc_types[boundary] == 1:
                    Df = k * area / L

                    rows.append(c1); cols.append(c1); data.append(Df)
                    b[c1] += Df * self.vars.bcs[boundary].val[i]

                # Neumann (currently commented in your code)
                elif self.params.bc_types[boundary] == 2:
                    b[c1] -= self.params.bcs[boundary].val[i] * area

        # build CSR matrix
        A = coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

        return A, b