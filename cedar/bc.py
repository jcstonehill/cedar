from abc import ABC, abstractmethod

import cedar


class BC(ABC):
    def __init__(self):
        self.boundary: str = None
        self.mesh: cedar.Mesh = None

    @abstractmethod
    def initialize(self, t_start: float):
        pass

    @abstractmethod
    def update(self, t: float):
        pass