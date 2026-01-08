#https://ntrs.nasa.gov/citations/20240008350

import numpy as np

import cedar

class Parahydrogen(cedar.base.Fluid):
    """
    Parahydrogen

    Reference:

    TODO Cite Parahydrogen Properties v05
    """

    name = "Parahydrogen"

    def __init__(self):
        raise Exception("Not implemented.")
    
    def rho_from_T_P(self, T: float, P: float) -> float:
        pass
    
    def cp_from_T_P(self, T: float, P: float) -> float:
        pass

    def mu_from_T_P(self, T: float, P: float) -> float:
        pass

    def k_from_T_P(self, T: float, P: float) -> float:
        pass

    def e_from_T_P(self, T: float, P: float) -> float:
        pass

    def T_from_e_P(self, e: float, P: float) -> float:
        pass