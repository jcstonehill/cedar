import matplotlib.pyplot as plt
import numpy as np
import os

import cedar


class DirichletSteady1D(cedar.Benchmark):
    """
    ``Heat Transfer`` benchmark.
    """

    is_quick: bool = True
    tol: float = 2

    @classmethod
    def compute(cls) -> dict:
        mesh = cedar.Mesh3D(os.path.dirname(os.path.abspath(__file__)) + "/data/dirichlet_steady_1d.vtkhdf")
        
        material = cedar.materials.ConstantMaterial(1, 1, 1)

        ht = cedar.HeatTransfer(mesh)
        ht.materials = {"volume" : material}
        ht.bc = {
            "left" : cedar.FixedTemperatureBC(0),
            "right" : cedar.FixedTemperatureBC(1)
        }
        ht.T.ic = 0.5

        problem = cedar.Problem([ht])
        problem.output = None

        problem.solve()

        x0 = mesh.cell_centers[:, 0]
        T0 = ht.T.as_continuous_cell_value()
        sort_indices = np.argsort(x0)

        x = x0[sort_indices]
        T = T0[sort_indices]

        T_ref = np.copy(x)

        return {"T" : T}, {"T" : T_ref}, {"x_pos" : x}

    @classmethod
    def plot(cls, show = True, save = False):
        out, ref, meta = cls.compute()

        plt.plot(meta["x_pos"], out["T"], linewidth = 1)
        
        plt.plot(meta["x_pos"], ref["T"], ":", linewidth = 1)

        plt.grid(True)
        plt.title("Heat Transfer Benchmark: Dirichlet Steady 1D")
        plt.xlabel("x")
        plt.ylabel("Temperature [K]")
        plt.legend(["Cedar", "Analytical"])

        if save:
            plt.savefig(f"{str(cls.__name__).lower()}_T.png")

        if show:
            plt.show()

        plt.clf()