from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt


class Material(ABC):
    """
    Abstract base class for solid material properties.

    This interface defines temperature-dependent properties commonly
    used in heat transfer and neutron transport calculations.
    """
    
    full_name = "Material"

    T_min_k = 0
    T_max_k = 10000

    T_min_cp = 0
    T_max_cp = 10000

    def cp(self, T: np.ndarray) -> np.ndarray:
        """
        Compute specific heat capacity as a function of temperature.

        Temperature is clipped between T_min_cp and T_max_cp.

        :param T: Cell temperature array.
        :return: Specific heat capacity.
        """
        T = np.clip(T, self.T_min_cp, self.T_max_cp)
        return self._cp(T)

    def k(self, T: np.ndarray) -> np.ndarray:
        """
        Compute thermal conductivity as a function of temperature.

        Temperature is clipped between T_min_k and T_max_k.

        :param T: Cell temperature array.
        :return: Thermal conductivity.
        """
        T = np.clip(T, self.T_min_k, self.T_max_k)
        return self._k(T)

    @abstractmethod
    def _cp(self, T: np.ndarray) -> np.ndarray:
        """
        Compute specific heat capacity as a function of temperature.

        Temperature is clipped between T_min_cp and T_max_cp.

        :param T: Cell temperature array.
        :return
        """
        pass

    @abstractmethod
    def _k(self, T: np.ndarray) -> np.ndarray:
        """
        Compute thermal conductivity as a function of temperature.

        Temperature is clipped between T_min_k and T_max_k.

        :param T: Cell temperature array.
        :return: Thermal conductivity.
        """
        pass

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

    def plot_k(self, path: str, show = True, save = False):
        """
        Plot thermal conductivity over a standard temperature range.
        """
        T = np.linspace(0, 3500, 3500)

        y_lim = 1.1*max([self.k(self.T_min_k), self.k(3500)])
        plt.plot(T, self.k(T))
        plt.axvspan(self.T_min_k, self.T_max_k, alpha = 0.2)
        plt.xlim([0, 3500])
        plt.ylim([0, y_lim])
        plt.grid(True)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Thermal Conductivity [W/m-K]")
        plt.legend(["k Value", "Valid Temperature Range"])
        plt.title(self.__class__.__name__ + ": k")

        if show:
            plt.show()

        if save:
            plt.savefig(path + "/" + self.__class__.__name__.lower() + "_k.png")

        plt.clf()

    def plot_cp(self, path: str, show = True, save = False):
        """
        Plot specific heat capacity over a standard temperature range.
        """
        T = np.linspace(0, 3500, 3500)

        y_lim = 1.1*max([self.cp(self.T_min_cp), self.cp(3500)])
        plt.plot(T, self.cp(T))
        plt.axvspan(self.T_min_cp, self.T_max_cp, alpha = 0.2)
        plt.xlim([0, 3500])
        plt.ylim([0, y_lim])
        plt.grid(True)
        plt.xlabel("Temperature [K]")
        plt.ylabel("Specific Heat Capacity [J/kg-K]")
        plt.legend(["cp Value", "Valid Temperature Range"])
        plt.title(self.__class__.__name__ + ": cp")
        plt.savefig("material_plots/" + self.__class__.__name__.lower() + "_cp.png")

        if show:
            plt.show()

        if save:
            plt.savefig(path + "/" + self.__class__.__name__.lower() + "_cp.png")

        plt.clf()