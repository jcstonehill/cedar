import numpy as np
import matplotlib.pyplot as plt

import cedar


class NAFEMSTransient1D(cedar.base.Benchmark):
    """
    ``Thermal`` benchmark.
    """
    is_quick = True
    tol = 1

    @classmethod
    def plot(cls, show = True, save = False):
        data, ref_data, meta = cls._run()

        plt.plot(meta["x_pos"], meta["vals"], linewidth = 5)
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

    @classmethod
    def _run(cls):
        mesh = cedar.Mesh3D("cedar/benchmarks/thermal/nafems_transient_1d/nafems_transient_1d.msh")

        material = cedar.materials.ConstantMaterial(7200, 35, 440.5)

        thermal = cedar.models.Thermal("thermal", mesh)
        thermal.set_material("volume", material)
        thermal.set_bc("left", "dirichlet", 0)
        thermal.set_bc("right", "dirichlet", 0)
        thermal.vars.T.set_initial(0)
        thermal.vars.T.set(0)
        
        def right_bc_func(t):
            return 100 * np.sin(np.pi*t/40)
        
        right_bc = cedar.models.ScalarFromFunc("Right BC", "K", right_bc_func)

        problem = cedar.Problem("NAFEMS_1D_Transient", create_outputs = False)
        problem.add_model(right_bc)
        problem.add_model(thermal)

        problem.couple(right_bc.vars.scalar, thermal.vars.bcs["right"], cedar.adapters.NearestValue)

        problem.solve(dt = 1, t_end = 32)

        x_pos = np.array([pos[0] for pos in mesh.cell_centers])

        combined = zip(x_pos, thermal.vars.T.val)
        x_pos, vals = zip(*sorted(combined))

        data = {
            "point" : [np.interp(0.08, x_pos, vals)]
        }

        ref_data = {
            "point" : [36.6]
        }

        meta = {
            "x_pos" : x_pos,
            "vals" : vals
        }

        return data, ref_data, meta