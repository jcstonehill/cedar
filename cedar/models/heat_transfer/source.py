from abc import ABC, abstractmethod
import numpy as np
from typing import Callable, Iterable

import cedar


class HeatTransferSource(cedar.Term, ABC):
    """
    Abstract base class for :class:`cedar.HeatTransfer` heat sources.

    """

    @abstractmethod
    def _ht_contribution(self, cell_vols: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def _initialize(self) -> np.ndarray:
        pass


class HeatTransferTotalSource(cedar.Term):
    """
    Total heat generation source for :class:`cedar.HeatTransfer`.

    Distributes a total heating in [W] over the mesh using an optional spatial
    shape function to yield a heating source term.

    Parameters
    ----------
    Qgen : float | Callable
        Heat generation in [W].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face).
            - ex. Qgen = 1000
        - A function f(t) evaluated at simulation time.
            - ex. Qgen = lambda t: t + 1000

    shape : Callable
        Relative heat fraction shape. This is a function f(x, y, z) evaluated at
        cell center positions. Once evaluated, this is normalized such the sum =
        1 to ensure conservation of energy.

        ex. shape = lambda x, y, z: x + y + z

    """

    def __init__(
        self,
        Qgen: float | Callable = None,
        shape: Callable = None,
    ):
        self.Qgen: float | Callable = Qgen
        self.shape: Callable = shape

    def _ht_contribution(self, cell_vols: np.ndarray) -> np.ndarray:
        return self.values["Qgen"]

    def _initialize(self):
        # We have to do some manual manipulation since this has separation of
        # variables. The internal term._interpret_user_value() function can only
        # handle functions that take 4 arguments: x, y, z, t. We evaluate the
        # shape, normalize it, and then create a lambda function that multiplies
        # it by Qgen or Qgen(t)
        if self.shape is not None:
            x, y, z = self.pts[:, 0], self.pts[:, 1], self.pts[:, 2]
            shape = self.shape(x, y, z)

        else:
            shape = np.full(self.N, 1 / self.N, dtype=np.float64)

        shape = shape * self.geom_w
        shape = shape / np.sum(shape)

        if callable(self.Qgen):
            Qgen_func = lambda x, y, z, t: self.Qgen(t) * shape

        else:
            Qgen_func = lambda x, y, z, t: self.Qgen * shape

        self._add_value("Qgen", Qgen_func)


class HeatTransferVolumetricSource(cedar.Term):
    """
    Volumetric heat source for :class:`cedar.HeatTransfer`.

    Adds a volumetric heating in [W/m^3] as a heat source term.

    Parameters
    ----------
    volQgen : float | Iterable | Callable
        Heat generation in [W/m^3].

        Can be specified by one of the following:

        - A scalar value (same value enforced to each face). ex. volQgen = 1000
        - An iterable with size = number of cells (1 value per cell in order of cell_i). ex. volQgen = [1000, 1100, 1200, 1300, 1400] # where region has 5 cells.
        - A function f(x, y, z, t) evaluated at cell center positions and simulation time. ex. volQgen = lambda x, y, z, t: x + y + z + t + 1000

    Attributes
    ----------
    volQgen : float | Iterable | Callable
        Heat generation in [W/m^3].

    """

    def __init__(
        self,
        volQgen: float | Iterable | Callable = None,
    ):
        self.volQgen: float | Iterable | Callable = volQgen

    def _ht_contribution(self, cell_vols: np.ndarray) -> np.ndarray:
        return cell_vols * self.values["volQgen"]

    def _initialize(self):
        self._add_value("volQgen", self.volQgen)
