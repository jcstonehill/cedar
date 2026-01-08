import numpy as np

import cedar


class ZrC(cedar.base.Material):
    """
    Zirconium Carbide

    Reference:
    SNP-HDBK-0008 SNP Material Property Handbook
    https://ntrs.nasa.gov/citations/20240004217
    """

    def rho_rt(self):
        return 6730
    
    def k(self, T):
        # if T < 100 or T > 2650:
        #     raise Exception("Violated temperature limits: " + str(T))
           
        A0, A1, A2 = 23.76, 8.9, -0.7014

        return A0 + A1*(T/1000) + A2*(T/1000)**2
    
    def cp(self, T):
        T = np.asarray(T, dtype=float)
        T_k = T / 1000.0

        cp = np.empty_like(T_k)

        mask = T <= 294

        # Low-temperature branch
        N, A0, A1, A2 = 2.509, 52730, -4.986, 83.5
        cp[mask] = (
            A0 * T_k[mask]**N
            / (1 + A1*T_k[mask] + A2*T_k[mask]**2)
        )

        # High-temperature branch
        B0, B1, B2, B_2 = 488.9, -21.33, 29.64, -10.88
        cp[~mask] = (
            B0
            + B1*T_k[~mask]
            + B2*T_k[~mask]**2
            + B_2 / (T_k[~mask]**2)
        )

        return cp