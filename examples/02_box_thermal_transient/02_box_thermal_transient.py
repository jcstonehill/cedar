import cedar

# https://www.robertribando.com/modules/pipeflow/

# Instantiate Problem Object
problem = cedar.Problem("02_box_thermal_transient")

# Build Mesh for Thermal Model
mesh = cedar.Mesh3D("02_box_thermal_transient.msh")

# Define Material Properties
ZrC = cedar.materials.ZrC()

# Instantiate Model
thermal = cedar.models.Thermal("thermal", mesh, {"volume" : ZrC})

# Set Heat Source
thermal.set_Qdot(5e4)

# Set Boundary Conditions (Default is Adiabatic)
thermal.set_bc("left", "dirichlet", 300)
thermal.set_bc("right", "dirichlet", 500)

# Set Initial Conditions
thermal.vars.T.set_initial(500)

# Add Problem to Model and Solve
problem.add_model(thermal)
problem.solve(dt = 100, t_end = 10000)