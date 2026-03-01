import numpy as np

import cedar


class ConstantMaterial(cedar.Material):
    """
    Constant Material
    """

    full_name = "Constant Material"

    T_min_k = 0
    T_max_k = 10000

    T_min_cp = 0
    T_max_cp = 10000

    def __init__(self, rho_rt: float = 1, k: float = 1, cp: float = 1):
        self._rho_rt: float = rho_rt
        self._k_val: float = k
        self._cp_val: float = cp

    def rho_rt(self) -> float:
        return self._rho_rt
    
    def _k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        return np.full_like(T, self._k_val, dtype = np.float64)
    
    def _cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        return np.full_like(T, self._cp_val, dtype = np.float64)