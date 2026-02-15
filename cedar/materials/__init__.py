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

            T = np.linspace(material.T_min, material.T_max, 100)

            plt.plot(T, material.k(T))
            plt.grid(True)
            plt.xlabel("Temperature [K]")
            plt.ylabel("Thermal Conductivity [W/m-K]")
            plt.title(Material.__name__ + ": k")
            plt.savefig(path + "/" + Material.__name__.lower() + "_k.png")
            plt.clf()

            print("Created " + Material.__name__.lower() + "_k.png")

            plt.plot(T, material.cp(T))
            plt.grid(True)
            plt.xlabel("Temperature [K]")
            plt.ylabel("Specific Heat Capacity [J/kg-K]")
            plt.title(Material.__name__ + ": cp")
            plt.savefig(path + "/" + Material.__name__.lower() + "_cp.png")
            plt.clf()

            print("Created " + Material.__name__.lower() + "_cp.png")

        except:
            print("Couldn't plot " + Material.__name__)
