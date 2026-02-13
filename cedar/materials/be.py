import numpy as np

import cedar


class Be(cedar.Material):
    """
    Beryllium

    Reference:
        SNP-HDBK-0008 SNP Material Property Handbook
        https://ntrs.nasa.gov/citations/20240004217
    """

    T_min = 200
    T_max = 1560

    def rho_rt(self) -> float:
        return 1848
    
    def k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T / 1000.0
        T_k2 = T_k*T_k
        
        A0, A1, A2, AN, N = 182.6, -118.9, 20.47, 1.192, -2.825
        return A0 + A1*T_k + A2*T_k2 + AN*T_k**N
    
    def cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T / 1000.0
        T_k2 = T_k*T_k

        cp = np.zeros_like(T)

        mask = T <= 293

        A0, A1, A2 = 14230., -5.341, 14.39
        cp[mask] = (A0*T_k2[mask]) / (1 + A1*T_k[mask] + A2*T_k2[mask])

        # Coefficients multiplied by 1000 to convert from [J/g-K] to [J/kg-K]
        B0, B1, B2, B_2 = 2401., 337.5, 280.4, -63.88
        cp[~mask] = B0 + B1*T_k[~mask] + B2*T_k2[~mask] + B_2/T_k2[~mask]

        return cp        