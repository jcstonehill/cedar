import numpy as np

import cedar


class UO2(cedar.Material):
    """
    Uranium Dioxide

    Reference:
        J. V. Miller, “Estimating Thermal Conductivity of CERMET Fuel Materials
        for Nuclear Reactor Application,” Lewis Research Center, Report NASA TND-3898,
        Cleveland, OH, Apr. 1967.
        https://ntrs.nasa.gov/citations/19670013537

        S. G. Popov, et. al., "Thermophysical Properties of MOX and UO2 Fuel
        Including the Effects of Radiation," Oak Ridge National Laboratory, Report
        ORNL/TM-2000/351, Oak Ridge, TN, Nov. 2000.

    Notes:
        The temperature where the thermal conductivity correlation transitions
        from one curve fit to the other is shifted from 1800 K to 1620.54 K to
        yield a smooth transition.
    """

    full_name = "Uranium Dioxide"

    T_min_k = 500
    T_max_k = 3000
    
    T_min_cp = 298.15
    T_max_cp = 3120

    def __init__(self, P: float = 0):
        self.P = P

    def rho_rt(self) -> float:
        return (1-self.P)*10970.0
    
    def _k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T / 1000.0
        T_k2 = T_k*T_k
        k = np.zeros_like(T)

        mask = T <= 1620.54

        A0, A1, A2 = 10.41, -9.44, 2.52
        k[mask] = A0 + A1*T_k[mask] + A2*T_k2[mask]

        A0 = 1.73
        k[~mask] = A0

        return k
    
    def _cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T[T<=5] = 5
        
        C1, C2, C3, theta, Ea = 302.27, 8.463e-3, 8.741e7, 548.68, 18531.7
        theta_T = theta/T

        return C1*(theta_T)**2 * np.exp(theta_T)/((np.exp(theta_T)-1)**2) + 2*C2*T + C3*Ea*np.exp(-Ea/T)/(T**2)