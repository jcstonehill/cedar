from abc import ABC, abstractmethod
import numpy as np

import cedar


class Benchmark(ABC):
    """
    Abstract base class for validation benchmarks.

    A ``Benchmark`` defines a reproducible comparison between computed data and
    reference data, producing a single scalar error metric. Benchmarks also
    support plotting for diagnostic or reporting purposes.

    Attributes
    ----------
    is_quick : bool
        If this benchmark is considered to have low computation time.
    tol : float
        The acceptable MAPE for this benchmark to be considered passed.

    """

    is_quick: bool = True
    tol: float = None

    @classmethod
    @abstractmethod
    def compute(cls) -> tuple[dict, dict, dict]:
        """
        Build and solve the benchmark.

        Returns
        -------
        dict[str, np.ndarray]
            Calculated values from Cedar.
        dict[str, np.ndarray]
            Pre-defined reference values.
        dict[str, np.ndarray]
            All other values that might be useful. An example could be an array
            of positions assocaited with the values that can be used for
            plotting.

        """

        pass

    @classmethod
    @abstractmethod
    def plot(cls, show: bool = True, save: bool = False):
        """
        Get data from cls.compute() and then create relevant plots.

        Parameters
        ----------
        show : bool
            Show plot in new window.
        save : bool
            Save plots as png files.

        """

        pass

    @classmethod
    def print(cls):
        """
        Run the benchmark and print a formatted error summary.

        """

        error = cls.run()
        print(f"{cls.__name__} Error = {error:.1f} [%]")

    @classmethod
    def run(cls) -> float:
        """
        Execute the benchmark and compute the aggregate error.

        Returns
        -------
        float
            Mean Absolute Percentage Error between computed and reference data.

        """

        out, ref, meta = cls.compute()

        combined_out = []
        combined_ref = []

        for var in out:
            combined_out.extend(out[var])
            combined_ref.extend(ref[var])

        return cedar.functions.MAPE(combined_out, combined_ref)
