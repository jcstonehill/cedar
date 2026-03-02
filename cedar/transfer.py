from abc import ABC, abstractmethod
import numpy as np

import cedar


class FieldView:
    def __init__(self, field: "Field", key: str):
        self.field: "Field" = field
        self.key: str = key

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

    def as_continuous_cell_value(self) -> np.ndarray:
        value = np.zeros(self.mesh.N_cells)

        for region in self.mesh.regions:
            value[self.mesh.region_i[region]] = self.values[region]

        return value

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

    def ic_as_continuous_cell_value(self) -> np.ndarray:
        ic = np.zeros(self.mesh.N_cells)

        for region in self.mesh.regions:
            ic[self.mesh.region_i[region]] = self.ic[region]

        return ic

    def initialize(self, t_start: float):
        # Convert user inpute to dict[domain] = np.ndarray
        self.ic = self._interpret_user_ic(t_start)

        # Now, copy IC to be initial guess of values
        for key in self.ic:
            self.values[key] = np.copy(self.ic[key])

    def step(self, t: float):
        for key in self.ic:
            self.ic[key] = np.copy(self.values[key])

    def _interpret_user_ic(self, t_start: float) -> np.ndarray:
        ic = {}

        if self.is_on_regions:
            for region in self.mesh.regions:
                N = self.mesh.region_N[region]

                ic_input = self.ic[region] if isinstance(self.ic, dict) else self.ic

                # No IC Provided
                if ic_input is None:
                    ic[region] = np.zeros(N)

                # IC is func(x, y, z)
                elif callable(ic_input):
                    # Extract indices for the whole region at once
                    cell_indices = self.mesh.region_i[region]
                    
                    # Get coordinates for all cells using advanced indexing: Shape (N, 3)
                    coords = self.mesh.cell_centers[cell_indices]
 
                    ic[region] = np.array(ic_input(coords[:, 0], coords[:, 1], coords[:, 2], t_start), dtype=np.float64)
                
                elif np.array(ic_input).size == 1:
                    ic[region] = np.full(N, ic_input, dtype=np.float64)

                else:
                    ic[region] = np.copy(ic_input)

        if self.is_on_boundaries:
            for boundary in self.mesh.boundaries:
                N = self.mesh.boundary_N[boundary]

                ic_input = self.ic[boundary] if isinstance(self.ic, dict) else self.ic

                # No IC Provided
                if ic_input is None:
                    ic[boundary] = np.zeros(N)

                # IC is func(x, y, z)
                elif callable(ic_input):
                    # Extract indices for the whole region at once
                    face_indicies = self.mesh.boundary_i[boundary]
                    
                    # Get coordinates for all cells using advanced indexing: Shape (N, 3)
                    coords = self.mesh.face_centers[face_indicies]
 
                    ic[boundary] = np.array(ic_input(coords[:, 0], coords[:, 1], coords[:, 2], t_start), dtype=np.float64)
                
                elif np.array(ic_input).size == 1:
                    ic[boundary] = np.full(N, ic_input, dtype=np.float64)

                else:
                    ic[boundary] = np.copy(ic_input)

        return ic
    
    def __getitem__(self, key: str) -> FieldView:
        return FieldView(self, key)

class Transfer:
    def __init__(self, src: cedar.Field | cedar.FieldView, mapping: np.ndarray):
        self.src: cedar.Field | cedar.FieldView = src
        self.mapping: np.ndarray = mapping

    def get(self) -> np.ndarray:
        if isinstance(self.src, cedar.Field):
            return self.mapping @ self.src.as_continuous_cell_value()
        
        else:
            return self.mapping @ self.src.values()