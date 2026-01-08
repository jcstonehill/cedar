import numpy as np
from dataclasses import dataclass
from typing import Callable

import cedar

@dataclass
class Params:
    func: Callable

@dataclass
class Vars:
    scalar: cedar.ScalarVar

class ScalarFromFunc(cedar.base.Model):
    """
    Calculate scalar value from user func.
    """
    @property
    def mesh(self) -> None:
        return self._mesh
    
    @property
    def params(self) -> Params:
        return self._params
    
    @property
    def vars(self) -> Vars:
        return self._vars

    def __init__(self, name, units, func = None):
        self.name = name
        self._mesh = None

        self._params = Params(func = func)

        self._vars = Vars(
            scalar = self.add_var(cedar.ScalarVar(self.name, "scalar", 0, units))
        )
    
    def check(self):
        pass

    def setup(self):
        pass

    def iterate(self, t0 = None, dt = None, res_reduc = 1e-2):
        if t0 is None:
            self.vars.scalar.set(self.params.func())

        else:
            self.vars.scalar.set(self.params.func(t0+dt))

        return 0