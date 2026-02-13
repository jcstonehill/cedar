import numpy as np

import cedar


class ZrC(cedar.Material):
    """
    Zirconium Carbide

    Reference:
    SNP-HDBK-0008 SNP Material Property Handbook
    https://ntrs.nasa.gov/citations/20240004217
    """

    T_min = 100
    T_max = 2650

    def rho_rt(self) -> float:
        return 6730.0
    
    def k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T/1000.0
        T_k2 = T_k*T_k

        A0, A1, A2 = 23.76, 8.9, -7.014
        return A0 + A1*T_k + A2*T_k2
    
    def cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T/1000.0
        T_k2 = T_k*T_k
        cp = np.zeros_like(T_k)

        mask = T <= 293

        # Coefficients multiplied by 1000 to convert from [J/g-K] to [J/kg-K]
        N, A0, A1, A2 = 2.509, 52730, -4.986, 83.5
        cp[mask] = (
            A0 * T_k[mask]**N
            / (1 + A1*T_k[mask] + A2*T_k2[mask])
        )

        B0, B1, B2, B_2 = 488.9, -21.33, 29.64, -10.88
        cp[~mask] = (
            B0
            + B1*T_k[~mask]
            + B2*T_k2[~mask]
            + B_2 / (T_k2[~mask])
        )