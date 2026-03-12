import cedar

# Create Mesh for Heat Transfer Model
mesh = cedar.Mesh3D(
    "examples/02_box_transient_heat_transfer/box_transient_heat_transfer.vtkhdf"
)

# Instantiate Model
ht = cedar.HeatTransfer(mesh)

ht.materials = {"volume": cedar.materials.ZrC()}

# Set Heat Source
ht.sources = cedar.HeatTransferTotalSource(5e4)

# Set Boundary Conditions (Default is Adiabatic)
ht.bcs = {"left": cedar.FixedTemperatureBC(300), "right": cedar.FixedTemperatureBC(500)}

# Set Initial Conditions
ht.T.ic = 400

# Create Problem
problem = cedar.Problem()
problem.dt = 100
problem.t_end = 10000
problem.add_model(ht)

problem.solve()
