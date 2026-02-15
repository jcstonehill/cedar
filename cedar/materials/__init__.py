import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from cedar.materials.be import Be
from cedar.materials.beo import BeO
from cedar.materials.constant import ConstantMaterial
from cedar.materials.u10mo import U10Mo
from cedar.materials.uc_zrc_c import UC_ZrC_C
from cedar.materials.un import UN
from cedar.materials.uo2 import UO2
from cedar.materials.zrc import ZrC
from cedar.materials.zrc_c import ZrC_C

all_materials = [
    Be,
    BeO,
    ConstantMaterial,
    U10Mo,
    UC_ZrC_C,
    UN,
    UO2,
    ZrC,
    ZrC_C
]

def plot_all(path = "material_plots"):
    if os.path.exists(path):
        shutil.rmtree(path)

    os.makedirs(path)

    for Material in all_materials:
        try:
            material = Material()

            T = np.linspace(0, 3500, 3500)

            y_lim = 1.1*max([material.k(material.T_min_k), material.k(3500)])
            plt.plot(T, material.k(T))
            plt.axvspan(material.T_min_k, material.T_max_k, alpha = 0.2)
            plt.xlim([0, 3500])
            plt.ylim([0, y_lim])
            plt.grid(True)
            plt.xlabel("Temperature [K]")
            plt.ylabel("Thermal Conductivity [W/m-K]")
            plt.legend(["k Value", "Valid Temperature Range"])
            plt.title(Material.__name__ + ": k")
            plt.savefig(path + "/" + Material.__name__.lower() + "_k.png")
            plt.clf()

            print("Created " + Material.__name__.lower() + "_k.png")

            y_lim = 1.1*max([material.cp(material.T_min_cp), material.cp(3500)])
            plt.plot(T, material.cp(T))
            plt.axvspan(material.T_min_cp, material.T_max_cp, alpha = 0.2)
            plt.xlim([0, 3500])
            plt.ylim([0, y_lim])
            plt.grid(True)
            plt.xlabel("Temperature [K]")
            plt.ylabel("Specific Heat Capacity [J/kg-K]")
            plt.legend(["cp Value", "Valid Temperature Range"])
            plt.title(Material.__name__ + ": cp")
            plt.savefig(path + "/" + Material.__name__.lower() + "_cp.png")
            plt.clf()

            print("Created " + Material.__name__.lower() + "_cp.png")

        except:
            print("Couldn't plot " + Material.__name__)
