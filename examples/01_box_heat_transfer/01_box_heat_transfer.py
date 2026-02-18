import cedar

# Create Mesh for Heat Transfer Model
mesh = cedar.Mesh3D("examples/01_box_heat_transfer/box_heat_transfer.vtkhdf")

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
ht.T.ic = 400

# Create Problem
problem = cedar.Problem()
problem.add_model(ht)

problem.solve()