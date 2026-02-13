from abc import ABC, abstractmethod
import numpy as np

import cedar

class Model(ABC):
    # @property
    # def mesh(self) -> cedar.Mesh:
    #     return self._mesh

    # @mesh.setter
    # def mesh(self, new_mesh: cedar.Mesh):
    #     if hasattr(self, "_mesh"):
    #         raise Exception("Can't change model's mesh after instantiation.")

    #     self._mesh = new_mesh

    def __init__(self):
        self.mesh: cedar.Mesh = None
        self.bc: dict[str, cedar.BC] = {}
        self.fields: dict[str, cedar.Field] = {}
        self.source: cedar.Source = None
        self.output_fields: list[str] = []
        
    def add_field(self, name: str, requires_ic = False, is_on_regions: bool = True, is_on_boundaries: bool = False) -> cedar.Field:
        if not hasattr(self, "fields"):
            self.fields: dict[str, cedar.Field] = {}
            
        if not hasattr(self, "output_fields"):
            self.output_fields: list[str] = []

        new_field = cedar.Field(name, requires_ic, is_on_regions, is_on_boundaries)
        self.fields[name] = new_field

        self.output_fields.append(name)

        return new_field

    # @abstractmethod
    # def check(self):
    #     pass

    @abstractmethod
    def initialize(self, t_start: float):
        pass

    @abstractmethod
    def iterate(self, res_reduc: float, dt: float) -> float:
        pass