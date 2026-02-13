import numpy as np

import cedar


class BeO(cedar.Material):
    """
    Beryllium Oxide

    Reference:
        SNP-HDBK-0008 SNP Material Property Handbook
        https://ntrs.nasa.gov/citations/20240004217
    """

    T_min = 200
    T_max = 2301

    def rho_rt(self) -> float:
        return 3010
    
    def k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T / 1000.0
        T_k2 = T_k*T_k
        
        A0, A1, A_0, A_1 = 57.1, 0.3584, 0.0476, 0.2395
        return (A0 + A1*T_k) / (A_0 + A_1*T_k + T_k2)
    
    def cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T / 1000.0
        T_k2 = T_k*T_k

        cp = np.zeros_like(T)

        mask = T <= 293

        N, A0, A1, A2 = 2.715, 47440., -4.321, 22.82
        cp[mask] = (A0*T_k[mask]**N) / (1 + A1*T_k[mask] + A2*T_k2[mask])

        # Coefficients multiplied by 1000 to convert from [J/g-K] to [J/kg-K]
        B0, B1, B_2 = 1645, 347.3, -66.13
        cp[~mask] = B0 + B1*T_k[~mask] + B_2/T_k2[~mask]

        return cp        