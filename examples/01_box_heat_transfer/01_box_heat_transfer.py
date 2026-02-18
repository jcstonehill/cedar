import os

import cedar

# Create Mesh for Heat Transfer Model
mesh = cedar.Mesh3D("examples/01_box_thermal/box_heat_transfer.vtkhdf")

# Instantiate Model
ht = cedar.HeatTransfer(mesh)

ht.materials = {
    "volume" : cedar.materials.ZrC()
}

# Set Heat Source
ht.source = cedar.HeatSource(5e4)

# Set Boundary Conditions (Default is Adiabatic)
ht.bc = {
    "left" : cedar.FixedTemperatureBC(300),
    "right" : cedar.FixedTemperatureBC(500)
}

# Set Initial Conditions
ht.T.ic = 500

# Create Problem
problem = cedar.Problem()
problem.models.append(ht)

problem.solve()

os.system(f"cp output.vtkhdf /mnt/c/Users/jacob/Documents/output.vtkhdf")