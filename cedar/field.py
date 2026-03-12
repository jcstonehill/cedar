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


class Field:
    def __init__(
        self,
        name: str,
        model_name: str,
        default_val: float,
        is_on_regions: bool = True,
        is_on_boundaries: bool = False,
    ):
        self.name: str = name
        self.model_name: str = model_name
        self.default_val: float = default_val
        self.is_on_regions: bool = is_on_regions
        self.is_on_boundaries: bool = is_on_boundaries

        self.ic: np.ndarray = (
            None  # Could be a number, a dict of numbers, a function, or a dict of functions
        )
        self.mesh: cedar.Mesh = (
            None  # This is assigned externally while the problem is built
        )

        self._values: dict[str, np.ndarray] = {}
        self._ic: dict[str, np.ndarray] = {}

    def get(self, domain: str = None):
        if domain is None:
            values = np.zeros(self.mesh.N_cells, dtype=np.float64)

            for region in self.mesh.regions:
                values[self.mesh.region_i[region]] = self._values[region]

            return values

        return self._values[domain]

    def set_ic(self, ic):
        self.ic = ic

    def _get_ic(self, domain: str = None):
        if self._ic is None:
            cedar.Log.error(
                f'Field named "{self.model_name}.{self.name}" requires an initial condition.'
            )

        if domain is None:
            ic = np.zeros(self.mesh.N_cells, dtype=np.float64)

            for region in self.mesh.regions:
                ic[self.mesh.region_i[region]] = self._ic[region]

            return ic

        return self._ic[domain]

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

    def _initialize(self):
        # Convert user input to dict[domain] = np.ndarray
        if self.ic is not None:
            self._ic = self._interpret_ic_from_user(self.ic)

            for key in self._ic:
                self._values[key] = np.copy(self._ic[key])

            if self.is_on_boundaries:
                cell_values = self.get()

                for boundary in self.mesh.boundaries:
                    boundary_i = self.mesh.boundary_i[boundary]
                    self._values[boundary] = np.copy(
                        cell_values[self.mesh.face_cells_i[boundary_i, 0]]
                    )

        else:
            if self.is_on_regions:
                for region in self.mesh.regions:
                    self._values[region] = np.full(
                        self.mesh.region_N[region], self.default_val, dtype=np.float64
                    )

            if self.is_on_boundaries:
                for boundary in self.mesh.boundaries:
                    self._values[boundary] = np.full(
                        self.mesh.boundary_N[boundary],
                        self.default_val,
                        dtype=np.float64,
                    )

    def _set(self, values: np.ndarray, domain: str = None):
        if domain is None:
            for region in self.mesh.regions:
                self._values[region][:] = values[self.mesh.region_i[region]]

        else:
            self._values[domain][:] = values

    def _step(self, t: float):
        if self._ic is not None:
            for key in self._ic:
                self._ic[key][:] = self._values[key]

    def _interpret_ic_from_user(self, user_ic) -> np.ndarray:
        ic = {}

        if not self.is_on_regions:
            cedar.Log.error(
                f"You can not supply an initial condition for a boundary-only field: {self.model_name}.{self.name}"
            )

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

                ic[region] = np.array(
                    ic_input(coords[:, 0], coords[:, 1], coords[:, 2]), dtype=np.float64
                )

            elif np.array(ic_input).size == 1:
                ic[region] = np.full(N, ic_input, dtype=np.float64)

            else:
                ic[region] = np.copy(ic_input)

        return ic

    def __getitem__(self, key: str) -> FieldView:
        return FieldView(self, key)
