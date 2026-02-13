from abc import ABC, abstractmethod
import numpy as np

import cedar


class Source(ABC):
    def __init__(self):
        self.mesh: cedar.Mesh = None

    @abstractmethod
    def initialize(self, t_start: float):
        pass

    @abstractmethod
    def iterate(self, T: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        pass

    @abstractmethod
    def step(self, t: float):
        pass