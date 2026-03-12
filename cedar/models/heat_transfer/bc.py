from abc import ABC, abstractmethod
from typing import Callable, Iterable
import numpy as np

import cedar


class HeatTransferBC(cedar.Term, ABC):
    """
    Abstract base class for :class:`cedar.HeatTransfer` boundary conditions.
    """

    @abstractmethod
    def _boundary_J(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        """
        Compute the heat flux through boundary faces.

        """

        pass

    @abstractmethod
    def _boundary_T(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        """
        Compute the temperature at boundary faces.

        """

        pass

    @abstractmethod
    def _ht_contribution(
        self, T_face: np.ndarray, k: np.ndarray, d: np.ndarray, A: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Return contributions for A, b matrices of :class:`cedar.HeatTransfer`.

        """
        pass


class AdiabaticBC(HeatTransferBC):
    """
    Adiabatic (zero-flux) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces zero normal heat flux across the boundary.
    """

    def _boundary_J(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        return np.zeros_like(T_cell)

    def _boundary_T(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        return np.copy(T_cell)

    def _initialize(self):
        pass

    def _ht_contribution(
        self, T_face: np.ndarray, k: np.ndarray, d: np.ndarray, A: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        LHS = np.zeros(self.N, dtype=np.float64)
        RHS = np.zeros(self.N, dtype=np.float64)

        return LHS, RHS


class KnownTemperatureBC(HeatTransferBC):
    """
    Known-temperature (Dirichlet) boundary condition for
    :class:`cedar.HeatTransfer`.

    Enforces a prescribed boundary temperature.

    Parameters
    ----------
    T : float or Iterable or f(x, y, z, t) or cedar.Field or cedar.FieldView
        Known temperature to apply to boundary in [K].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face).
            - ex. T = 300
        - An iterable with size = number of faces (1 value per face in order of face_i).
            - ex. T = [300, 400, 500, 600, 700] # where boundary has 5 faces.
        - A function f(x, y, z, t) evaluated at face center positions and simulation time.
            - ex. T = lambda x, y, z, t: x + y + z + t + 300
        - A cedar.Field/cedar.FieldView which will be used to create model coupling using a Transfer with nearest value mapping.
            - ex. T = flow.T_wall

    Attributes
    ----------
    T : float or Iterable or f(x, y, z, t) or cedar.Field or cedar.FieldView
        Known temperature to apply to boundary in [K].

    """

    def __init__(self, T: float | Iterable | Callable | cedar.Field | cedar.FieldView):
        self.T: float | Iterable | Callable | cedar.Field | cedar.FieldView = T

    def _boundary_J(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        T = self.values["T"]

        return k * (T_cell - T) / d

    def _boundary_T(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        return self.values["T"]

    def _initialize(self):
        T = self.T

        if self._is_field(T):
            T = cedar.Transfer.with_nearest_value_mapping(T, self.pts)

        self._add_value("T", T)

    def _ht_contribution(
        self, T_face: np.ndarray, k: np.ndarray, d: np.ndarray, A: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        # T is provided as an argument but we use our bc value here instead
        T = self.values["T"]

        LHS = k * A / d
        RHS = T * k * A / d

        return LHS, RHS


class KnownFluxBC(HeatTransferBC):
    """
    Known-flux (Neumann) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces a prescribed boundary flux.

    Parameters
    ----------
    J : float or Iterable or f(x, y, z, t)
        Known flux to apply to boundary in [W/m^2].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face).
            - ex. J = 1
        - An iterable with size = number of faces (1 value per face in order of face_i).
            - ex. J = [1, 2, 3, 4, 5] # where boundary has 5 faces.
        - A function f(x, y, z, t) evaluated at face center positions and simulation time.
            - ex. J = lambda x, y, z, t: x + y + z + t

    Attributes
    ----------
    J : float or Iterable or f(x, y, z, t)
        Known flux to apply to boundary in [W/m^2].

    """

    def __init__(self, J: float | Callable):
        self.J: float | Callable = J

    def _boundary_J(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        return self.values["J"]

    def _boundary_T(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        J = self.values["J"]

        return T_cell - d * J / k

    def _initialize(self):
        self._add_value("J", self.J)

    def _ht_contribution(
        self, T_face: np.ndarray, k: np.ndarray, d: np.ndarray, A: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        J = self.values["J"]

        LHS = np.zeros(self.N, dtype=np.float64)
        RHS = -J * A

        return LHS, RHS


class ConvectiveBC(HeatTransferBC):
    """
    Convective (Robin) boundary condition for :class:`cedar.HeatTransfer`.

    Enforces a boundary flux based on convection to a fluid, J = htc*(T_face - T_flow).

    Parameters
    ----------
    T_flow : float or Iterable or f(x, y, z, t) or cedar.Field or cedar.FieldView
        Flow temperature at the boundary face in [K].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. T_flow = 300
        - An iterable with size = number of faces (1 value per face in order of face_i). ex. T_flow = [300, 400, 500, 600, 700] # where boundary has 5 faces.
        - A function f(x, y, z, t) evaluated at face center positions and simulation time. ex. T_flow = lambda x, y, z, t: x + y + z + t + 300
        - A cedar.Field/cedar.FieldView which will be used to create model coupling using a Transfer with nearest value mapping. ex. T_flow = flow.T

    htc : float or Iterable or f(x, y, z, t) or cedar.Field or cedar.FieldView
        Heat transfer coefficient in [W/m^2-K].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. htc = 1000.
        - An iterable with size = number of faces (1 value per face in order of face_i). ex. htc = [1000, 1100, 1200, 1300, 1400] # where boundary has 5 faces.
        - A function f(x, y, z, t) evaluated at face center positions and simulation time. ex. htc = lambda x, y, z, t: x + y + z + t + 1000.
        - A cedar.Field/cedar.FieldView which will be used to create model coupling using a Transfer with nearest value mapping. ex. htc = flow.htc

    Attributes
    ----------
    T_flow : float or Iterable or f(x, y, z, t) or cedar.Field or cedar.FieldView
        Flow temperature at the boundary face in [K].
    htc : float or Iterable or f(x, y, z, t) or cedar.Field or cedar.FieldView
        Heat transfer coefficient in [W/m^2-K].

    """

    def __init__(
        self,
        T_flow: float | Callable | cedar.FieldView,
        htc: float | Callable | cedar.FieldView,
    ):
        self.T_flow: float | Callable | cedar.FieldView = T_flow
        self.htc: float | Callable | cedar.FieldView = htc

    def _boundary_J(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        T_flow = self.values["T_flow"]
        htc = self.values["htc"]

        return k * htc * (T_cell - T_flow) / (d * htc + k)

    def _boundary_T(
        self, T_cell: np.ndarray, k: np.ndarray, d: np.ndarray
    ) -> np.ndarray:
        T_flow = self.values["T_flow"]
        htc = self.values["htc"]

        return (k * T_cell + htc * d * T_flow) / (htc * d + k)

    def _initialize(self):
        T_flow = self.T_flow
        htc = self.htc

        if self._is_field(T_flow):
            T_flow = cedar.Transfer.with_nearest_value_mapping(T_flow, self.pts)

        if self._is_field(self.htc):
            htc = cedar.Transfer.with_nearest_value_mapping(htc, self.pts)

        self._add_value("T_flow", T_flow)
        self._add_value("htc", htc)

    def _ht_contribution(
        self, T_face: np.ndarray, k: np.ndarray, d: np.ndarray, A: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        T_flow = self.values["T_flow"]
        htc = self.values["htc"]

        C = A * htc * k / (d * htc + k)

        LHS = C
        RHS = C * T_flow

        return LHS, RHS
