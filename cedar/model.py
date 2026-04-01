from abc import ABC, abstractmethod

import cedar


class Model(ABC):

    @property
    def mesh(self) -> cedar.Mesh:
        return self._mesh

    @mesh.setter
    def mesh(self, new_mesh: cedar.Mesh):
        self._error("Unable to set new mesh after model instantiation.")

    model_name = "Model"

    def __init__(self, mesh: cedar.Mesh):
        self._mesh: cedar.Mesh = mesh

        self.scalar_fields: dict[str, cedar.ScalarField] = {}
        self.boundary_face_fields: dict[str, cedar.BoundaryFaceField] = {}
        self.cell_fields: dict[str, cedar.CellField] = {}
        self.all_fields: dict[str, cedar.Field] = {}

        self.output_fields: list[str] = []

        self.bcs: dict[str, cedar.Term] = {}
        self.sources: dict[str, cedar.Term] = {}

    def _add_boundary_face_field(
        self,
        name: str,
        default_val: float,
    ) -> cedar.BoundaryFaceField:
        new_field = cedar.BoundaryFaceField(name, self.model_name, default_val)

        self.boundary_face_fields[name] = new_field
        self.all_fields[name] = new_field
        self.output_fields.append(name)

        return new_field

    def _add_cell_field(
        self,
        name: str,
        default_val: float,
    ) -> cedar.CellField:
        new_field = cedar.CellField(name, self.model_name, default_val)

        self.cell_fields[name] = new_field
        self.all_fields[name] = new_field
        self.output_fields.append(name)

        return new_field

    def _add_scalar_field(
        self,
        name: str,
        default_val: float,
    ) -> cedar.ScalarField:
        new_field = cedar.ScalarField(name, self.model_name, default_val)

        self.scalar_fields[name] = new_field
        self.all_fields[name] = new_field
        self.output_fields.append(name)

        return new_field

    def _error(self, message: str):
        cedar.Log.error(f"{self.model_name}: " + message)

    @abstractmethod
    def _build(self):
        pass

    @abstractmethod
    def _initialize(self):
        pass

    @abstractmethod
    def _iterate(self, res_reduc: float, dt: float) -> float:
        pass
