import numpy as np

import cedar


class YH188(cedar.Material):
    """
    Yttrium Hydride (YH1.88)

    Reference:
        X. Hu, et. al., “Handbook on the Material Properties of Yttrium Hydride
        for High-Temperature Moderator Applications,” Oak Ridge National
        Laboratory, Report ORNL/TM-2021/2052, Oak Ridge, TN, Jun. 2021.
    """

    full_name = "Yttrium Hydride (YH1.88)"

    T_min_k = 300
    T_max_k = 823

    T_min_cp = 300
    T_max_cp = 823

    def rho_rt(self) -> float:
        return 4260.584
    
    def a(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T[T<300] = 300
        T[T>823] = 823

        return 1/(160*T-24600.)

    def _k(self, T: np.ndarray) -> np.ndarray:
        return self.rho_rt()*self.a(T)*self.cp(T)
    
    def _cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T[T<300] = 300
        T[T>823] = 823

        T2 = T*T
        cp = np.zeros_like(T)

        mask = T <= 648

        # Coefficients multiplied by 1000 to convert from [J/g-K] to [J/kg-K]
        a, b, c = 53.5, 1.2, -2.26e-4
        cp[mask] = a+b*T[mask]+c*T2[mask]

        # Coefficients multiplied by 1000 to convert from [J/g-K] to [J/kg-K]
        a, b, c, d = 2.23e13, -4.09e-2, 350, 9.88e-4
        cp[~mask] = a*np.exp(b*T[~mask]) + c*np.exp(d*T[~mask])

        return cp        