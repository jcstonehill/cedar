import numpy as np
from abc import ABC, abstractmethod
from collections.abc import Iterable


class Fluid(ABC):
    """
    Abstract base class defining a thermophysical fluid properties.

    This interface specifies the minimum set of thermodynamic and transport
    property relations required for a fluid to be used in Cedar. All properties
    are expressed as functions of temperature and pressure unless otherwise
    noted.

    Notes
    -----
    - Implementations may use analytical correlations, tabulated data, equations
      of state, or external libraries.
    - All methods are assumed to be valid over a specific range of temperatures
      and pressures, which should be documented by each concrete implementation.
    - All methods can be evaluated with a single point or iterables.
    - All units are SI unless otherwise specified.
    """

    @abstractmethod
    def rho_from_T_P(self, T: np.ndarray, P: np.ndarray) -> np.ndarray:
        """
        Compute fluid density from temperature and pressure.

        Parameters
        ----------
        T : ndarray
            Temperature.
        P : ndarray
            Pressure.

        Returns
        -------
        ndarray:
            Fluid density.
        """
        pass

    @abstractmethod
    def cp_from_T_P(self, T: np.ndarray, P: np.ndarray) -> np.ndarray:
        """
        Compute specific heat capacity at constant pressure.

        Parameters
        ----------
        T : ndarray
            Temperature.
        P : ndarray
            Pressure.

        Returns
        -------
        ndarray
            Specific heat capacity at constant pressure.
        """
        pass

    @abstractmethod
    def mu_from_T_P(self, T: np.ndarray, P: np.ndarray) -> np.ndarray:
        """
        Compute dynamic viscosity from temperature and pressure.

        Parameters
        ----------
        T : ndarray
            Temperature.
        P : ndarray
            Pressure.

        Returns
        -------
        ndarray
            Dynamic viscosity.
        """
        pass

    @abstractmethod
    def k_from_T_P(self, T: np.ndarray, P: np.ndarray) -> np.ndarray:
        """
        Compute thermal conductivity from temperature and pressure.

        Parameters
        ----------
        T : ndarray
            Temperature.
        P : ndarray
            Pressure.

        Returns
        -------
        ndarray
            Thermal conductivity.
        """
        pass

    @abstractmethod
    def e_from_T_P(self, T: np.ndarray, P: np.ndarray) -> np.ndarray:
        """
        Compute specific internal energy from temperature and pressure.

        Parameters
        ----------
        T : ndarray
            Temperature.
        P : ndarray
            Pressure.

        Returns
        -------
        ndarray
            Specific internal energy.
        """
        pass

    @abstractmethod
    def T_from_e_P(self, e: np.ndarray, P: np.ndarray) -> np.ndarray:
        """
        Compute temperature from specific internal energy and pressure.

        Parameters
        ----------
        e : ndarray
            Specific internal energy.
        P : ndarray
            Pressure.

        Returns
        -------
        ndarray
            Temperature.
        """
        pass