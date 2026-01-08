import numpy as np
import matplotlib.pyplot as plt

import cedar


class AdiabaticTransient0D(cedar.base.Benchmark):
    """
    ``Thermal`` benchmark.
    """
    is_quick = True
    tol = 1

    @classmethod
    def plot(cls, show = True, save = False):
        data, ref_data, meta = cls._run()

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

    @classmethod
    def _run(cls):
        mesh = cedar.Mesh3D("cedar/benchmarks/thermal/adiabatic_transient_0d/adiabatic_transient_0d.msh")

        material = cedar.materials.ConstantMaterial(rho_rt = 1, k = 1, cp = 1)

        thermal = cedar.models.Thermal("thermal", mesh)
        thermal.set_material("volume", material)
        thermal.set_Qdot(mesh.vol)
        thermal.vars.T.set_initial(0)
        thermal.vars.T.set(0)
    
        problem = cedar.Problem("Adiabatic Transient 0D", create_outputs = False)
        problem.add_model(thermal)

        T_vals = np.zeros(26)

        for i in range(25):
            problem.solve(dt = 1, t0 = i, t_end = i+1)
            T_vals[i+1] = np.average(thermal.vars.T.val)

        data = {
            "T" : T_vals
        }

        ref_data = {
            "T" : np.linspace(0., 25., 26)
        }

        meta = {
            "t_array" : np.linspace(0., 25., 26)
        }
        
        return data, ref_data, meta