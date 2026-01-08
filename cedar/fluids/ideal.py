import numpy as np

import cedar


class IdealGas(cedar.base.Fluid):
    """
    Ideal Gas
    """
    name = "IdealGas"

    def __init__(self, k: float, cp: float, mu: float, molar_mass: float):
        self.k = k
        self.cp = cp
        self.mu = mu
        self.molar_mass = molar_mass

        self.R = 8314.46261815324 / molar_mass
        self.cv = self.cp - self.R

    def rho_from_T_P(self, T: float, P: float) -> float:
        return P/(self.R*T)
    
    def cp_from_T_P(self, T: float, P: float) -> float:
        return np.full_like(T, self.cp)

    def mu_from_T_P(self, T: float, P: float) -> float:
        return np.full_like(T, self.mu)

    def k_from_T_P(self, T: float, P: float) -> float:
        return np.full_like(T, self.k)

    def e_from_T_P(self, T: float, P: float) -> float:
        return self.cv*T

    def T_from_e_P(self, e: float, P: float) -> float:
        return e/self.cv