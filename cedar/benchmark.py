from abc import ABC, abstractmethod
import numpy as np

import cedar


class Benchmark(ABC):
    """
    Abstract base class for validation benchmarks.

    A ``Benchmark`` defines a reproducible comparison between computed data
    and reference data, producing a single scalar error metric. Benchmarks
    also support plotting for diagnostic or reporting purposes.
    """

    is_quick: bool = True
    tol: float = None

    @classmethod
    @abstractmethod
    def compute(cls) -> tuple[dict, dict]:
        pass

    @classmethod
    @abstractmethod
    def plot(cls, show: bool = True, save: bool = False):
        pass

    @classmethod
    def print(cls):
        """
        Run the benchmark and print a formatted error summary.
        """
        error = cls.run()
        print(f"{cls.__name__} Error = {error:.1f} [%]")

    @classmethod
    def run(cls) -> bool:
        """
        Execute the benchmark and compute the aggregate error.

        Returns
        -------
        bool
            Benchmark passed

        """
        out, ref, meta = cls.compute()

        combined_out = []
        combined_ref = []

        for var in out:
            combined_out.extend(out[var])
            combined_ref.extend(ref[var])
        
        return cedar.functions.MAPE(combined_out, combined_ref)