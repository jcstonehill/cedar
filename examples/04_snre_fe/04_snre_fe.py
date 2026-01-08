import cedar
import numpy as np

# Instantiate Problem
problem = cedar.Problem("04_snre_fe")

# THERMAL MODEL
# =============================================================================
# Qdot for 1/6th of a FE
fe_Qdot = 0.94809*366.97e6/(564*6)

# Calculate heat flux to tie tube element
Qdot_to_tie_tube = 0.061*fe_Qdot
flux_to_tie_tube = Qdot_to_tie_tube / 0.00978911

# Instantiate Model
thermal_mesh = cedar.Mesh3D("04_snre_fe.msh")
thermal = cedar.models.Thermal("thermal", thermal_mesh, {"fuel": cedar.materials.UC_ZrC_C()})

# Set the global heating value for the thermal model
thermal.set_Qdot(fe_Qdot)

# Set a cosine relative heating fraction shape.
thermal.vars.Qdot_shapes["fuel"].set([np.sin(np.pi*center[2]/0.89) for center in thermal_mesh.cell_centers])

# Change the BC for wall to the tie tube to heat flux instead of adiabatic.
thermal.set_bc("wall_to_tietube", "neumann", flux_to_tie_tube)

# Calculate ICs
T_ic = [370 + (3000-370)*center[2]/0.89 for center in thermal_mesh.cell_centers]
T_wall_ic = [370 + (3000-370)*center[2]/0.89 for center in thermal_mesh.face_centers[thermal_mesh.boundaries["fc_wall"]]]

# Set the BC at fc_wall to Dirichlet for coupling.
thermal.set_bc("fc_wall", "dirichlet", T_wall_ic)

# Set better initial guess to speed up convergence.
thermal.vars.T.set_initial(T_ic)

# FLOW MODEL
# =============================================================================
flow_mesh = cedar.Mesh1D(200, 0, 0.89, "fc", "fc_wall")
fluid = cedar.fluids.IdealGas(k = 0.983, cp = 17500, mu = 0.0000334346923392326, molar_mass = 2.01588)
fc_or = 0.23654/200
flow = cedar.models.Flow("flow", flow_mesh, fluid, Dh = 2*fc_or, A = (np.pi*fc_or**2)*19/6, P_wall = np.pi*2*fc_or*19/6, eps = 25e-6)
 
# Flow BCs
flow.vars.inlet.set([370.1, 3.1e6, 8.33/(564*6)])
flow.vars.outlet.set([370.1, 3.1e6, 8.33/(564*6)])

# Flow ICs
flow.vars.P.set_initial(3.1e6)
flow.vars.T.set_initial(370)
flow.vars.u.set_initial(10)

# COUPLING
# =============================================================================
problem.couple(thermal.vars.boundary_Qdots["fc_wall"], flow.vars.Qdot, cedar.adapters.Summation)
problem.couple(thermal.vars.boundary_Qdots["fc_wall"], flow.vars.Qdot_shape, cedar.adapters.Summation)
problem.couple(flow.vars.T_wall, thermal.vars.bcs["fc_wall"], cedar.adapters.NearestValue)

# SOLVE
# =============================================================================
problem.add_model(thermal)
problem.add_model(flow)

problem.solve()