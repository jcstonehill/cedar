import numpy as np

import cedar


class FieldView:
    def __init__(self, field: "Field", key: str):
        self.field: "Field" = field
        self.key: str = key

        self.N: int = self.field.values[key].size

    def values(self) -> np.ndarray:
        return self.field.values[self.key]
    
    @property
    def mesh(self) -> cedar.Mesh:
        return self.field.mesh

class Field:
    def __init__(
            self,
            name: str,
            requires_ic: bool = False,
            is_on_regions: bool = True,
            is_on_boundaries: bool = False
    ):
        self.name: str = name
        self.requires_ic: bool = requires_ic
        self.is_on_regions: bool = is_on_regions
        self.is_on_boundaries: bool = is_on_boundaries

        self.ic: np.ndarray  = None         # Could be a number, a dict of numbers, a function, or a dict of functions
        self.mesh: cedar.Mesh = None        # This is assigned externally while the problem is built

        self.values: dict[str, np.ndarray] = {}

    def get_ic(self, domain: str = None):
        if domain is None:
            ic = np.zeros(self.mesh.N_cells, dtype = np.float64)

            for region in self.mesh.regions:
                ic[self.mesh.region_i[region]] = self.ic[region]

            return ic
        
        return self.ic[domain]

    def get(self, domain: str = None):
        if domain is None:
            values = np.zeros(self.mesh.N_cells, dtype = np.float64)

            for region in self.mesh.regions:
                values[self.mesh.region_i[region]] = self.values[region]

            return values
        
        return self.values[domain]
    
    # def check(self):
    #     pass
        # if self.requires_ic and self.ic is None:
        #     raise Exception(f"Field requires IC: {self.name}")

        # if not self.requires_ic 

        # if isinstance(self.ic, dict):
        #     if self.is_on_regions:
        #         for region in self.mesh.regions:
        #             if region not in self.ic:
        #                 raise Exception(f"No IC provided for region: {region}")

    def initialize(self, t_start: float):
        # Convert user input to dict[domain] = np.ndarray
        self.ic = self._interpret_ic_from_user(self.ic)

        # Now, copy IC to be initial guess of values
        for key in self.ic:
            self.values[key] = np.copy(self.ic[key])

        if self.is_on_boundaries:
            for boundary in self.mesh.boundaries:
                self.values[boundary] = np.zeros(self.mesh.boundary_N[boundary], dtype = np.float64)

    def set(self, values: np.ndarray, domain: str = None):
        if domain is None:
            for region in self.mesh.regions:
                self.values[region][:] = values[self.mesh.region_i[region]]

        else:
            self.values[domain][:] = values

    def step(self, t: float):
        for key in self.ic:
            self.ic[key][:] = self.values[key]

    def _interpret_ic_from_user(self, user_ic) -> np.ndarray:
        ic = {}

        if self.is_on_regions:
            for region in self.mesh.regions:
                N = self.mesh.region_N[region]

                ic_input = user_ic[region] if isinstance(user_ic, dict) else user_ic

                # No IC Provided
                if ic_input is None:
                    ic[region] = np.zeros(N)

                # IC is func(x, y, z)
                elif callable(ic_input):
                    # Extract indices for the whole region at once
                    cell_indices = self.mesh.region_i[region]
                    
                    # Get coordinates for all cells using advanced indexing: Shape (N, 3)
                    coords = self.mesh.cell_centers[cell_indices]
 
                    ic[region] = np.array(ic_input(coords[:, 0], coords[:, 1], coords[:, 2]), dtype=np.float64)
                
                elif np.array(ic_input).size == 1:
                    ic[region] = np.full(N, ic_input, dtype=np.float64)

                else:
                    ic[region] = np.copy(ic_input)

        return ic
    
    def __getitem__(self, key: str) -> FieldView:
        return FieldView(self, key)