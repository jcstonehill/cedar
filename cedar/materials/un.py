import numpy as np

import cedar


class UN(cedar.Material):
    """
    Uranium Nitride

    Reference:
        SNP-HDBK-0008 SNP Material Property Handbook
        https://ntrs.nasa.gov/citations/20240004217
    """

    full_name = "Uranium Nitride"

    T_min_k = 10
    T_max_k = 2500

    T_min_cp = 5
    T_max_cp = 3000

    def __init__(self, P: float = 0):
        self.P = P

    def rho_rt(self) -> float:
        return (1-self.P)*14330.0
    
    def _k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T / 1000.0
        k = np.zeros_like(T)

        mask = T <= 2025

        A0 = 22.55
        N = 0.3845
        k[mask] = A0 * T_k[mask]**N

        A0 = 29.58
        k[~mask] = A0

        return k
    
    def _cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T / 1000.0
        T_k2 = T_k*T_k
        T_k3 = T_k2*T_k
        cp = np.zeros_like(T)

        mask = T <= 293

        N, A0, A1, A2, A3 = 2.168, 2.784e5, 9.191e1, 7.568e2, 4.086e2
        cp[mask] = (A0*T_k[mask]**N) / (1 + A1*T_k[mask] + A2*T_k2[mask] + A3*T_k3[mask])

        # Coefficients multiplied by 1000 to convert from [J/g-K] to [J/kg-K]
        B0, B1, B2, B_2 = 214.1, 9.746, 17.18, -2.607
        cp[~mask] = B0 + B1*T_k[~mask] + B2*T_k2[~mask] + B_2/T_k2[~mask]

        return cp        