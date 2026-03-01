import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from cedar.materials.be import Be
from cedar.materials.beo import BeO
from cedar.materials.constant import ConstantMaterial
from cedar.materials.g348 import G348
from cedar.materials.u10mo import U10Mo
from cedar.materials.uc_zrc_c import UC_ZrC_C
from cedar.materials.un import UN
from cedar.materials.uo2 import UO2
from cedar.materials.yh188 import YH188
from cedar.materials.zrc import ZrC
from cedar.materials.zrc_c import ZrC_C

all_materials = [
    Be,
    BeO,
    ConstantMaterial,
    G348,
    U10Mo,
    UC_ZrC_C,
    UN,
    UO2,
    YH188,
    ZrC,
    ZrC_C
]

def plot_all(path = "material_plots"):
    if os.path.exists(path):
        shutil.rmtree(path)

    os.makedirs(path)

    for Material in all_materials:
        # try:
        #     material = Material()
        material = Material()
        material.plot_k(path, show = False, save = True)
        material.plot_cp(path, show = False, save = True)

        # except:
        #     print("Couldn't plot " + Material.__name__)
