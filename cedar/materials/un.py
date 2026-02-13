import numpy as np

import cedar


class UN(cedar.Material):
    """
    Uranium Nitride

    Reference:
    SNP-HDBK-0008 SNP Material Property Handbook
    https://ntrs.nasa.gov/citations/20240004217
    """

    def __init__(self, P: float = 0):
        self.P = P

    def rho_rt(self) -> float:
        return (1-self.P)*14330.0
    
    def k(self, T: np.ndarray) -> np.ndarray:
        # if T < 100 or T > 2650:
        #     raise Exception("Violated temperature limits: " + str(T))
           
        A0, A1, A2 = 23.76, 8.9, -0.7014

        return A0 + A1*(T/1000) + A2*(T/1000)**2
    
    def cp(self, T: np.ndarray) -> np.ndarray:
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

        raise Exception("Need to implement this like zrc_c")