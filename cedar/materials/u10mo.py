import numpy as np

import cedar


class U10Mo(cedar.Material):
    """
    Uranium Molybdenum Alloy (10% weight Mo)

    Reference:
        D. E. Burkes, et. al., "Thermophysical Properties of U-10Mo Alloy," Idaho
        National Laboratory, Report INL/EXT-10-19373, Idaho Falls, ID, Nov. 2010.
    """

    T_min = 373.15
    T_max = 1073.15

    def __init__(self, P: float = 0):
        self.P = P

    def rho_rt(self) -> float:
        return (1-self.P)*17150.0
    
    def k(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)

        return 0.612435+0.0351*T
    
    def cp(self, T: np.ndarray) -> np.ndarray:
        T = np.array(T, dtype = np.float64)
        T2 = T*T

        C0, C1, C2 = 137, 0.0512, 0.0000199
        return C0 + C1*T + C2*T2

        