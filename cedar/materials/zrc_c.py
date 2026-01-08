import cedar
import numpy as np


class ZrC_C(cedar.base.Material):
    """
    Zirconium Carbide - Graphite Composite

    Reference:
    SNP-HDBK-0008 SNP Material Property Handbook
    https://ntrs.nasa.gov/citations/20240004217
    """

    def rho_rt(self):
        return 6700
    
    def k(self, T):
        T = np.asarray(T, dtype=float)
        T2 = T*T

        return np.where(
            T <= 1100,
            80.2218845 + (-0.067901)*T + 0.000022781*T2,
            24.012     + 0.008175*T   + 7.545e-8*T2
        )

    def cp(self, T):
        T = np.asarray(T, dtype=float)
        T2 = T*T
        T3 = T*T*T

        return 228.03 + 0.4422*T + -0.0001948*T2 + 3.318e-8*T3