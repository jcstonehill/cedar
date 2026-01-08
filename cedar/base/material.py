from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from collections.abc import Iterable


class Material(ABC):
    """
    Abstract base class for solid material properties.

    This interface defines temperature-dependent properties commonly
    used in heat transfer and neutron transport calculations.
    """

    _plot_T = np.linspace(0, 3000, 500)

    @abstractmethod
    def rho_rt(self) -> float:
        """
        Return material density at room temperature.

        Returns
        -------
        float
            Density at reference (room) temperature.
        """
        pass

    @abstractmethod
    def k(self, T: np.ndarray) -> np.ndarray:
        """
        Compute thermal conductivity as a function of temperature.

        Parameters
        ----------
        T : ndarray
            Temperature.

        Returns
        -------
        ndarray
            Thermal conductivity.
        """
        pass

    @abstractmethod
    def cp(self, T: np.ndarray) -> np.ndarray:
        """
        Compute specific heat capacity as a function of temperature.

        Parameters
        ----------
        T : ndarray
            Temperature.

        Returns
        -------
        ndarray
            Specific heat capacity.
        """
        pass

    def plot_k(self):
        """
        Plot thermal conductivity over a standard temperature range.
        """
        plt.plot(self._plot_T, self.k(self._plot_T))
        plt.grid(True)
        plt.xlabel("T [K]")
        plt.ylabel("k [W/m-K]")
        plt.show()
        plt.clf()

    def plot_cp(self):
        """
        Plot specific heat capacity over a standard temperature range.
        """
        plt.plot(self._plot_T, self.cp(self._plot_T))
        plt.grid(True)
        plt.xlabel("T [K]")
        plt.ylabel("Cp [J/kg-K]")
        plt.show()
        plt.clf()