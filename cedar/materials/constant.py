import numpy as np

import cedar


class ConstantMaterial(cedar.Material):
    """
    Constant Material
    """
    def __init__(self, rho_rt: float, k: float, cp: float):
        self._rho_rt: float = rho_rt
        self._k: float = k
        self._cp: float = cp

    def rho_rt(self) -> float:
        return self._rho_rt
    
    def k(self, T: np.ndarray) -> np.ndarray:
        return np.full_like(T, self._k)
    
    def cp(self, T: np.ndarray) -> np.ndarray:
        return np.full_like(T, self._cp)