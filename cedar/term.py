from abc import ABC, abstractmethod
import numpy as np
from typing import Callable

import cedar


class Term(ABC):
    @abstractmethod
    def _initialize(self, mesh: cedar.Mesh):
        pass

    def _add_value(self, name: str, user_value: np.ndarray | Callable | cedar.Transfer):
        self.values[name] = np.zeros(self.N, dtype=np.float64)
        self._values_update_funcs[name] = self._interpret_user_value(user_value)

    def _build(self, pts: np.ndarray, N: float, geom_w: np.ndarray):
        self.pts: np.ndarray = pts
        self.N: float = N
        self.geom_w: np.ndarray = geom_w / np.sum(geom_w)

        self.values: dict[str, np.ndarray] = {}
        self._values_update_funcs: dict[str, Callable] = {}

        self.transfers: list[cedar.Transfer] = []

    def _interpret_user_value(
        self, user_value: np.ndarray | Callable | cedar.Transfer
    ) -> Callable:
        """
        Convert user value input to a Callable f(t) for self.update()

        """
        if callable(user_value):
            x, y, z = self.pts[:, 0], self.pts[:, 1], self.pts[:, 2]
            return lambda t: user_value(x, y, z, t)

        elif isinstance(user_value, cedar.Transfer):
            self.transfers.append(user_value)
            return lambda t: user_value.get()

        # This can be either a number or an array since its assigned with
        # _value[:] = user_value
        return lambda t: user_value

    def _is_field(self, user_value) -> bool:
        return isinstance(user_value, cedar.Field) or isinstance(
            user_value, cedar.FieldView
        )

    def _update(self, t: float):
        for name in self.values:
            self.values[name][:] = self._values_update_funcs[name](t)
