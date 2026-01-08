from abc import ABC, abstractmethod
import numpy as np

class Benchmark(ABC):
    """
    Abstract base class for validation benchmarks.

    A ``Benchmark`` defines a reproducible comparison between computed data
    and reference data, producing a single scalar error metric. Benchmarks
    may optionally support plotting for diagnostic or reporting purposes.
    """

    is_quick = True
    tol = None

    @classmethod
    @abstractmethod
    def plot(cls, show: bool = True, save: bool = False):
        """
        Plot benchmark results.

        Parameters
        ----------
        show : bool, optional
            If True, display the plot interactively.
        save : bool, optional
            If True, save the plot to disk.
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
            Mean absolute percentage error (MAPE) in percent.

        Notes
        -----
        This method aggregates all benchmark data into a single error value
        using :meth:`MAPE`.
        """
        data, ref_data, _ = cls._run()

        combined_data = []
        combined_ref_data = []

        for key in data.keys():
            combined_data.extend(data[key])
            combined_ref_data.extend(ref_data[key])

        return cls.MAPE(combined_data, combined_ref_data)

    @classmethod
    def MAPE(cls, val, ref) -> float:
        """
        Compute the mean absolute percentage error (MAPE).

        Parameters
        ----------
        val : ndarray
            Computed or predicted values.
        ref : ndarray
            Reference values.

        Returns
        -------
        float
            Mean absolute percentage error in percent.

        Notes
        -----
        Zero values are replaced internally to avoid division by zero.
        """
        val = np.array(val)
        ref = np.array(ref)

        # Prevent division by zero
        ref[ref == 0] = 1e-64
        val[val == 0] = 1e-64

        return float(np.sum(np.abs((val - ref) / ref) * 100) / val.size)

    @classmethod
    @abstractmethod
    def _run(cls) -> tuple[dict, dict, dict]:
        """
        Execute the benchmark and return raw results.

        Returns
        -------
        data : dict
            Computed benchmark data, keyed by quantity name.
        ref_data : dict
            Reference benchmark data, keyed by quantity name.
        meta : dict
            Optional metadata associated with the benchmark run.
        """
        return {}, {}, {}