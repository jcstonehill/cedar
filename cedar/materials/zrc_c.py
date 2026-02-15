import cedar
import numpy as np


class ZrC_C(cedar.Material):
    """
    Zirconium Carbide - Graphite Composite

    Reference:
        SNP-HDBK-0008 SNP Material Property Handbook
        https://ntrs.nasa.gov/citations/20240004217

    Notes:
        The temperature where the thermal conductivity correlation transitions
        from one curve fit to the other is shifted from 1100 K to 1019.34 K to
        yield a smooth transition.
    """

    T_min_k = 200
    T_max_k = 3100

    T_min_cp = 300
    T_max_cp = 3200

    def rho_rt(self) -> float:
        return 6700.0
    
    def k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T/1000.0
        T_k2 = T_k * T_k
        k = np.zeros_like(T)

        mask = T <= 1019.34

        A0, A1, A2 = 77.967, -67.901, 22.781
        k[mask] = A0 + A1*T_k[mask] + A2*T_k2[mask]

        # A0 (below) was reduced from 24.012 to 21.757 to create a smooth curve.
        # It appears that this would still line up with the experimental data
        # used in SNP-HDBK-0008
        A0, A1, A2 = 24.012, 8.175, 0.07545
        k[~mask] = A0 + A1*T_k[~mask] + A2*T_k2[~mask]

        return k

    def cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T_k = T/1000.0
        T_k2 = T_k*T_k
        T_k3 = T_k2*T_k

        # Coefficients multiplied by 1000 to convert from [J/g-K] to [J/kg-K]
        A0, A1, A2, A3 = 228.03, 442.2, -194.8, 33.18
        return A0 + A1*T_k + A2*T_k2 + A3*T_k3