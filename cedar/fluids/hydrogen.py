from CoolProp.CoolProp import PropsSI
import numpy as np

import cedar

class Hydrogen(cedar.base.Fluid):
    """
    Hydrogen

    Reference:

    PropsSI
    """


    name = "Hydrogen"

    def rho_from_T_P(self, T: float, P: float) -> float:
        return PropsSI("D", "T", T, "P", P, self.name)
    
    def cp_from_T_P(self, T: float, P: float) -> float:
        return PropsSI("C", "T", T, "P", P, self.name)

    def mu_from_T_P(self, T: float, P: float) -> float:
        return PropsSI("V", "T", T, "P", P, self.name)

    def k_from_T_P(self, T: float, P: float) -> float:
        return PropsSI("L", "T", T, "P", P, self.name)

    def e_from_T_P(self, T: float, P: float) -> float:
        return PropsSI("U", "T", T, "P", P, self.name)

    def T_from_e_P(self, e: float, P: float) -> float:
        return PropsSI("T", "U", e, "P", P, self.name)
    
    def T_max(self):
        return PropsSI("TMAX", self.name)