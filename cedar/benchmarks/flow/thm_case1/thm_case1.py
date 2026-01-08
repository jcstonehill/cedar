import numpy as np
import matplotlib.pyplot as plt

import cedar


class THMCase1(cedar.base.Benchmark):
    """
    ``Flow`` benchmark.
    """
    is_quick = True
    tol = 1

    @classmethod
    def plot(cls, show = True, save = False):
        data, ref_data, meta = cls._run()

        x_pos = meta["x_pos"]

        y_labels = {
            "T" : "Temperature [K]",
            "P" : "Pressure [Pa]",
            "u" : "Velocity [m/s]",
            "rho" : "Density [kg/m3]",
            "e" : "Specific Internal Energy [J/kg]",
            "Re" : "Reynolds",
            "Pr" : "Prandtl"
        }

        for key in data.keys():
            plt.plot(x_pos, data[key], linewidth = 5)
            plt.plot(x_pos, ref_data[key], ":", linewidth = 3)
            plt.grid(True)
            plt.title("Flow Benchmark: THM Case 1")
            plt.xlabel("x [m]")
            plt.ylabel(y_labels[key])
            plt.legend(["Cedar", "THM"])
            plt.tight_layout()

            if save:
                plt.savefig(f"{str(cls.__name__).lower()}_{key}.png")

            if show:
                plt.show()

            plt.clf()

    @classmethod
    def _run(cls):
        problem = cedar.Problem("test_flow", create_outputs = False)

        mesh = cedar.Mesh1D(1000, 0, 1, "fc", "fc_wall")
        fluid = cedar.fluids.IdealGas(k = 0.7018, cp = 13803, mu = 2.8710e-5, molar_mass = 2.01588)
        flow = cedar.models.Flow("flow", mesh, fluid, Dh = 0.003, A = 0.00000706858, P_wall = 0.00942477796, eps = 200e-6)

        flow.vars.inlet.set([300, 5008026, 0.00077368421])
        flow.vars.outlet.set([300, 5008026, 0.00077368421])
        flow.vars.P.set_initial(5e6)

        flow.vars.Qdot.set(14137.16)

        problem.add_model(flow)
        problem.solve()

        x_pos = np.linspace(0, 1, 1000)

        data = {
            "T" : flow.vars.T.val,
            "P" : flow.vars.P.val,
            "u" : flow.vars.u.val,
            "rho" : flow.vars.rho.val,
            "e" : flow.vars.e.val,
            "Re" : flow.vars.Re.val,
            "Pr" : flow.vars.Pr.val
        }

        ref_ids = ["cp", "cv", "e", "k", "M", "mu", "P", "Pr", "Re", "rho", "rhoA", "rhoEA", "rhouA", "T", "v", "u"]
        ref_data = {}
        for ref_id in ref_ids:
            ref_data[ref_id] = []

        with open("cedar/benchmarks/flow/thm_case1/ref.csv") as file:
            lines = file.readlines()

            lines.pop(0)

            for line in lines:
                split = line.split(",")
                split.pop(0)

                for i in range(len(split)):
                    ref_data[ref_ids[i]].append(float(split[i]))

        meta = {
            "x_pos" : x_pos
        }

        return data, ref_data, meta