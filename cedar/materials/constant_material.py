import cedar


class ConstantMaterial(cedar.base.Material):
    """
    Generic Constant Material
    """
    def __init__(self, rho_rt, k, cp):
        self._rho_rt = rho_rt
        self._k = k
        self._cp = cp
        
    def rho_rt(self):
        return self._rho_rt
    
    def k(self, T):
        return self._k
    
    def cp(self, T):
        return self._cp