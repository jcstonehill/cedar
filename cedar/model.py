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
        self.fields: dict[str, cedar.Field] = {}
        self.output_fields: list[str] = []

        self.bcs: dict[str, cedar.Term] = {}
        self.sources: dict[str, cedar.Term] = None

    def _add_field(
        self,
        name: str,
        model_name: str,
        default_val: float,
        is_on_regions: bool = True,
        is_on_boundaries: bool = False,
    ) -> cedar.Field:
        if not hasattr(self, "fields"):
            self.fields: dict[str, cedar.Field] = {}

        if not hasattr(self, "output_fields"):
            self.output_fields: list[str] = []

        new_field = cedar.Field(
            name, model_name, default_val, is_on_regions, is_on_boundaries
        )
        self.fields[name] = new_field

        self.output_fields.append(name)

        return new_field

    # @abstractmethod
    # def check(self):
    #     pass

    @abstractmethod
    def _build(self):
        pass

    @abstractmethod
    def _initialize(self):
        pass

    @abstractmethod
    def _iterate(self, res_reduc: float, dt: float) -> float:
        pass
