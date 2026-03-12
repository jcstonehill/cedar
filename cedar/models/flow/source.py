from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import fsolve
from typing import Callable

import cedar


class FlowSource(cedar.Term, ABC):
    """
    Abstract base class for :class:`cedar.Flow` heat sources.

    """

    @abstractmethod
    def _initialize(self):
        pass

    @abstractmethod
    def _wall_heating_values(self, T, Re, Pr, k, Dh, A_wall, ds):
        pass


class FlowQdotSource(FlowSource):
    """
    Total heat source for :class:`cedar.Flow`.

    Distributes a total heating in [W] over the mesh using an optional spatial
    shape function to yield a heating source term.

    Parameters
    ----------
    Qdot : float | Callable
        Heat transfer in [W].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. Qdot = 1000
        - A function f(t) evaluated at simulation time. Qdot = lambda t: t + 1000

    shape : Callable
        Relative heat fraction shape. This is a function f(x, y, z) evaluated at
        cell center positions. Once evaluated, this is normalized such the sum =
        1 to ensure conservation of energy.

        ex. shape = lambda x, y, z: x + y + z

    """

    def __init__(
        self,
        Qdot: float | Callable = None,
        shape: Callable = None,
    ):
        self.Qdot: float | Callable = Qdot
        self.shape: Callable = shape

    def _initialize(self):
        # We have to do some manual manipulation since this has separation of
        # variables. The internal term._interpret_user_value() function can only
        # handle functions that take 4 arguments: x, y, z, t. We evaluate the
        # shape, normalize it, and then create a lambda function that multiplies
        # it by Qgen or Qgen(t)
        if self.shape is not None:
            x, y, z = self.pts[:, 0], self.pts[:, 1], self.pts[:, 2]
            shape = self.shape(x, y, z)
            shape = shape / np.sum(shape)

        else:
            shape = np.full(self.N, 1 / self.N, dtype=np.float64)

        if callable(self.Qdot):
            Qdot_func = lambda x, y, z, t: self.Qdot(t) * shape

        else:
            Qdot_func = lambda x, y, z, t: self.Qdot * shape

        self._add_value("Qdot", Qdot_func)

    def _wall_heating_values(self, T, Re, Pr, k, Dh, P_wall, ds):
        Qdot = self.values["Qdot"]

        s = np.arange(self.N) * ds + 0.5 * ds

        T_wall_func = (
            lambda T_wall: Qdot
            - cedar.functions.westinghouse_modified_mccarthy_wolf(
                Re, Pr, T_wall, T, s, Dh
            )
            * k
            / Dh
            * P_wall
            * ds
            * (T_wall - T)
        )

        T_wall = fsolve(T_wall_func, T)

        htc = cedar.functions.westinghouse_modified_mccarthy_wolf(
            Re, Pr, T_wall, T, s, Dh
        )
        Nu = htc * Dh / k

        return Qdot, T_wall, Nu, htc


class FlowHeatFluxSource(FlowSource):
    """
    Wall heat flux source term for :class:`cedar.HeatTransfer`.

    Adds a wall heat flux in [W/m^2] as a heat source term.

    Parameters
    ----------
    J : float | Iterable | Callable
        Wall heat flux in [W/m^2].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. J = 1000
        - An iterable with size = number of cells (1 value per cell in order of cell_i). ex. J = [1000, 1100, 1200, 1300, 1400] # where boundary has 5 faces.
        - A function f(x, y, z, t) evaluated at cell center positions and simulation time. ex. J = lambda x, y, z, t: x + y + z + t + 1000
        - A cedar.Field/cedar.FieldView which will be used to create model coupling using a Transfer with nearest value mapping. ex. J = heat_transfer.J["boundary"]

    Attributes
    ----------
    J : float | Iterable | Callable
        Wall heat flux in [W/m^2].

    """

    def __init__(self, J):
        self.J: float | Callable | cedar.FieldView = J

    def _initialize(self):
        if self._is_field(self.J):
            self.J = cedar.Transfer.with_nearest_value_mapping(self.J, self.pts)

        self._add_value("J", self.J)

    def _wall_heating_values(self, T, Re, Pr, k, Dh, P_wall, ds):
        J = self.values["J"]

        Qdot = P_wall * ds * J

        s = np.arange(self.N) * ds + 0.5 * ds

        T_wall_func = (
            lambda T_wall: Qdot
            - cedar.functions.westinghouse_modified_mccarthy_wolf(
                Re, Pr, T_wall, T, s, Dh
            )
            * k
            / Dh
            * P_wall
            * ds
            * (T_wall - T)
        )

        T_wall = fsolve(T_wall_func, T)

        htc = cedar.functions.westinghouse_modified_mccarthy_wolf(
            Re, Pr, T_wall, T, s, Dh
        )
        Nu = htc * Dh / k

        return Qdot, T_wall, Nu, htc


# class FlowWallTemperatureSource(FlowSource):
#     def __init__(self, T_wall: float | np.ndarray | Callable | cedar.FieldView):
#         self.T_wall = T_wall

#     def _initialize(self):
#         self._add_value("T_wall", self.T_wall)

#     def _wall_heating_values(self, T, Re, Pr, k, Dh, P_wall, ds):
#         T_wall = self.values["T_wall"]
#         s = np.arange(self.N) * ds + 0.5 * ds

#         Nu = cedar.functions.westinghouse_modified_mccarthy_wolf(
#             Re, Pr, self.T_wall, T, s, Dh
#         )
#         htc = Nu * k / Dh
#         Qdot = htc * P_wall * ds * (T_wall - T)

#         return Qdot, T_wall, Nu, htc
