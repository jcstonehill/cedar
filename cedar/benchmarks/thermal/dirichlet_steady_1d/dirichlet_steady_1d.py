import numpy as np
import matplotlib.pyplot as plt

import cedar


class DirichletSteady1D(cedar.base.Benchmark):
    """
    ``Thermal`` benchmark.
    """
    
    is_quick = True
    tol = 2.7

    @classmethod
    def plot(cls, show = True, save = False):
        data, ref_data, meta = cls._run()
        plt.plot(meta["x_pos"], data["T"], linewidth = 3)

        
        plt.plot(meta["x_pos"], ref_data["T"], ":", linewidth = 3)

        plt.grid(True)
        plt.title("Thermal Benchmark: Dirichlet Steady 1D")
        plt.xlabel("x")
        plt.ylabel("Temperature [K]")
        plt.legend(["Cedar", "Analytical"])

        if save:
            plt.savefig(f"{str(cls.__name__).lower()}_T.png")

        if show:
            plt.show()

        plt.clf()

    @classmethod
    def _run(cls):
        mesh = cedar.Mesh3D("cedar/benchmarks/thermal/dirichlet_steady_1d/dirichlet_steady_1d.msh")

        material = cedar.materials.ConstantMaterial(1, 1, 1)

        thermal = cedar.models.Thermal("thermal", mesh)
        thermal.set_material("volume", material)
        thermal.set_bc("left", "dirichlet", 0)
        thermal.set_bc("right", "dirichlet", 1)
        thermal.vars.T.set_initial(0.5)
        
        problem = cedar.Problem("DirichletSteady1D", create_outputs = False)
        problem.add_model(thermal)
        problem.solve()
        
        ref = np.zeros(mesh.N_cells)
        x_pos = np.zeros(mesh.N_cells)
        for i, (x, y, z) in enumerate(mesh.cell_centers):
            # Mesh is 0.1 [m] long and the values go from 0 to 1
            x_pos[i] = x
            ref[i] = x/0.1

        data = {
            "T" : thermal.vars.T.val
        }

        ref_data = {
            "T" : ref
        }

        meta = {
            "x_pos" : x_pos,
        }

        return data, ref_data, meta