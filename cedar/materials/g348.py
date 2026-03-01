import numpy as np

import cedar


class G348(cedar.Material):
    """
    G-348 Graphite

    Reference:
        D. M. McEligot, et. al., "Thermal Properties of G-348 Graphite,” Idaho
        National Laboratory, Report INL/EXT-16-38241, Idaho Falls, Idaho, May.
        2016.

        A. T. D. Butland, et. al., "The Specific Heat of Graphite: An Evaluation
        of Measurements," Journal of Nuclear Materials 49 (1973/74) 45-46, Jun.
        1973.
    """

    full_name = "G-348 Graphite"

    T_min_k = 303.15
    T_max_k = 1273.15

    T_min_cp = 200
    T_max_cp = 3500

    def rho_rt(self) -> float:
        return 1860.35
    
    def _k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T[T>self.T_max_k] = self.T_max_k
        T2 = T*T
        
        # Converted coeffs from C to K
        A0, A1, A2 = 166.111, -0.127717, 3.719e-5
        return A0 + A1*T + A2*T2
    
    def _cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T[T<self.T_min_cp] = self.T_min_cp
        T2 = T*T
        T3 = T2*T
        T4 = T3*T

        # Converted coefficients from cal/g-K to J/kg-K
        A0, A1 = 2269.748016, -0.010159981956
        A_1, A_2, A_3, A_4 = -3.77952e5, -1.819135e8, 6.6699492120e10, -6.015929184e12

        return A0 + A1*T + A_1/T + A_2/T2 + A_3/T3 + A_4/T4