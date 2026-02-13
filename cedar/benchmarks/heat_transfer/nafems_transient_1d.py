import matplotlib.pyplot as plt
import numpy as np
import os

import cedar


class NAFEMSTransient1D(cedar.Benchmark):
    """
    ``Heat Transfer`` benchmark.
    """

    is_quick: bool = True
    tol: float = 100

    @classmethod
    def compute(cls) -> dict:
        mesh = cedar.Mesh3D(os.path.dirname(os.path.abspath(__file__)) + "/data/nafems_transient_1d.vtkhdf")
        
        material = cedar.materials.ConstantMaterial(7200, 35, 440.5)

        bc_func = lambda x, y, z, t: 100*np.sin(np.pi*t/40)

        ht = cedar.HeatTransfer(mesh)
        ht.materials = {"volume" : material}
        ht.bc = {
            "left" : cedar.FixedTemperatureBC(0),
            "right" : cedar.FixedTemperatureBC(bc_func)
        }
        ht.T.ic = 0

        problem = cedar.Problem([ht])
        problem.dt = 1
        problem.t_end = 32
        problem.output = None

        problem.solve()

        x0 = mesh.cell_centers[:, 0]
        T0 = ht.T.as_continuous_cell_value()
        sort_indices = np.argsort(x0)

        x = x0[sort_indices]
        T = T0[sort_indices]

        return {"point" : [np.interp(0.08, x, T)]}, {"point" : [36.6]}, {"x_pos" : x, "T" : T}

    @classmethod
    def plot(cls, show = True, save = False):
        out, ref, meta = cls.compute()

        plt.plot(meta["x_pos"], meta["T"], linewidth = 5)
        plt.plot([0.08], [36.6], linestyle = "none", marker = "*", markersize = 14)
        plt.grid(True)
        plt.title("Thermal Benchmark: NAFEMS Transient 1D")
        plt.xlabel("x [m]")
        plt.ylabel("Temperature [C]")
        plt.legend(["Cedar", "NAFEMS"])

        if save:
            plt.savefig(f"{str(cls.__name__).lower()}_T.png")

        if show:
            plt.show()

        plt.clf()