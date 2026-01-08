import cedar

# Instantiate Problem Object
problem = cedar.Problem("03_heated_channel")

# Build 1D Mesh for Flow Model
mesh = cedar.Mesh1D(N_cells = 100, start_z = 0, end_z = 1, region = "fc", boundary = "fc_wall")

# Define Fluid Properties
fluid = cedar.fluids.IdealGas(k = 0.7018, cp = 13803, mu = 2.8710e-5, molar_mass = 2.01588)

# Instantiate Model
flow = cedar.models.Flow("flow", mesh, fluid, Dh = 0.003, A = 0.00000706858,
                         P_wall = 0.00942477796, eps = 200e-6)

# Set Heat Source
flow.vars.Qdot.set(14137.16)

# Set Flow Boundary Conditions
flow.vars.inlet.set([300, 5e6, 0.00077368421])
flow.vars.outlet.set([300, 5e6, 0.00077368421])

# Adjust initial conditions.
flow.vars.T.set_initial(300)
flow.vars.P.set_initial(5e6)
flow.vars.u.set_initial(1)

# Add Problem to Model and Solve
problem.add_model(flow)
problem.solve()