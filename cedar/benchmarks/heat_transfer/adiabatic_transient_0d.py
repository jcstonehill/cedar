import matplotlib.pyplot as plt
import numpy as np
import os

import cedar


class AdiabaticTransient0D(cedar.Benchmark):
    """
    ``Heat Transfer`` benchmark.
    """

    is_quick: bool = True
    tol: float = 1

    @classmethod
    def compute(cls) -> dict:
        mesh = cedar.Mesh3D(os.path.dirname(os.path.abspath(__file__)) + "/data/adiabatic_transient_0d.vtkhdf")
        
        material = cedar.materials.ConstantMaterial(1, 1, 1)

        ht = cedar.HeatTransfer(mesh)
        ht.materials = {"volume" : material}
        ht.T.ic = 0
        ht.source = cedar.HeatSource(volQdot=1)

        problem = cedar.Problem([ht])
        problem.output = None
        problem.t_end = 1

        T_vals = np.zeros(26)

        for i in range(25):
            problem.solve()

            T_vals[i+1] = np.average(ht.T.as_continuous_cell_value())

            problem.t_start += 1
            problem.t_end += 1

        return {"T" : T_vals}, {"T" : np.linspace(0., 25., 26)}, {"t_array" : np.linspace(0., 25., 26)}

    @classmethod
    def plot(cls, show = True, save = False):
        data, ref_data, meta = cls.compute()

        plt.plot(np.arange(26), data["T"], linewidth = 5)
        plt.plot(meta["t_array"], ref_data["T"], ":", linewidth = 3)
        plt.grid(True)
        plt.title("Thermal Benchmark: Adiabatic Transient 0D")
        plt.xlabel("t [s]")
        plt.ylabel("Temperature [K]")
        plt.legend(["Cedar", "Analytical"])

        if save:
            plt.savefig(f"{str(cls.__name__).lower()}_T.png")

        if show:
            plt.show()

        plt.clf()