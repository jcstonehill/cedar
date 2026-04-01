from abc import ABC, abstractmethod
import numpy as np

import cedar


class FieldView:
    """
    A view of a cedar.Field

    Contains a reference to a field and a domain. This can be treated as a field
    for coupling, but values from this will only come from the defined domain.

    Create a FieldView by using __getitem__ on a field:
    cedar.Field["domain"] -> FieldView

    Parameters
    ----------
    field : cedar.Field
        The source field.
    domain : str
        The domain from which to pull values.

    Attributes
    ----------
    field : cedar.Field
        The source field.
    domain : str
        The domain from which to pull values.
    mesh : cedar.Mesh
        Reference to source field's mesh.

    """

    def __init__(self, field: "Field", domain: str):
        self.field: "Field" = field
        self.domain: str = domain

    def get(self) -> np.ndarray:
        return self.field.get(self.domain)

    @property
    def mesh(self) -> cedar.Mesh:
        return self.field.mesh


class Field(ABC):
    def __init__(
        self,
        name: str,
        model_name: str,
        default_val: float,
    ):
        self.name: str = name
        self.model_name: str = model_name
        self.default_val: float = default_val

    @abstractmethod
    def get(self, domain: str) -> np.ndarray:
        pass

    @abstractmethod
    def _set(self, values: np.ndarray, domain: str = None):
        pass

    @abstractmethod
    def _initialize(self):
        pass

    def _build(self, mesh: cedar.Mesh):
        self.mesh: cedar.Mesh = mesh

    def __getitem__(self, key: str) -> FieldView:
        return FieldView(self, key)


class ScalarField(Field):
    def __init__(self, name: str, model_name: str, default_val: float):
        self.name: str = name
        self.model_name: str = model_name
        self.default_val: float = default_val

        self._values = None

    def get(self, domain: str) -> np.ndarray:
        return np.copy(self._values)

    def _get_ic(self, domain: str) -> np.ndarray:
        raise Exception("Initial condition not currently supported for scalar field.")

    def _set(self, values: np.ndarray, domain: str = None):
        self._values = values

    def _initialize(self):
        self._value = self.default_val


class BoundaryFaceField(Field):
    def __init__(self, name: str, model_name: str, default_val: float):
        self.name: str = name
        self.model_name: str = model_name
        self.default_val: float = default_val

        self._values: dict[str, np.ndarray] = {}

    def get(self, domain: str) -> np.ndarray:
        return np.copy(self._values[domain])

    def _initialize(self):
        self._values: dict[str, np.ndarray] = {}

        for boundary, N in self.mesh.boundary_N.items():
            self._values[boundary] = np.full(N, self.default_val, dtype=np.float64)

    def _set(self, values: np.ndarray, domain: str = None):
        self._values[domain][:] = values


class CellField(Field):
    def __init__(self, name: str, model_name: str, default_val: float):
        self.name: str = name
        self.model_name: str = model_name
        self.default_val: float = default_val

        # Could be a number, a dict of numbers, a function, or a dict of functions
        self.ic: np.ndarray = None

        self._values: np.ndarray = None
        self._ic: np.ndarray = None

    def get(self, domain: str = None) -> np.ndarray:
        if domain is None:
            return np.copy(self._values)

        return np.copy(self._values[self.mesh.region_i[domain]])

    def set_ic(self, ic):
        self.ic = ic

    def _get_ic(self, region: str = None):
        if self._ic is None:
            cedar.Log.error(
                f'Field named "{self.model_name}.{self.name}" requires an initial condition.'
            )

        if region is None:
            return self._ic

        return self._ic[self.mesh.region_i[region]]

    def _initialize(self):
        if self.ic is not None:
            self._ic = self._interpret_ic_from_user(self.ic)
            self._values = np.copy(self._ic)

        else:
            self._values = np.full(
                self.mesh.N_cells, self.default_val, dtype=np.float64
            )

    def _set(self, values: np.ndarray, domain: str = None):
        if domain is None:
            self._values[:] = values

        else:
            self._values[self.mesh.region_i[domain]] = values

    def _step(self):
        if self._ic is not None:
            self._ic[:] = self._values[:]

    def _interpret_ic_from_user(self, user_ic) -> np.ndarray:
        ic = np.zeros(self.mesh.N_cells, dtype=np.float64)

        for region in self.mesh.regions:
            region_i = self.mesh.region_i[region]

            ic_input = user_ic[region] if isinstance(user_ic, dict) else user_ic

            # IC is func(x, y, z)
            if callable(ic_input):
                region_i = self.mesh.region_i[region]

                x = self.mesh.cell_centers[region_i, 0]
                y = self.mesh.cell_centers[region_i, 1]
                z = self.mesh.cell_centers[region_i, 2]

                ic[region_i] = ic_input(x, y, z)

            # IC is scalar or iterable
            else:
                ic[region_i] = ic_input

        return ic
