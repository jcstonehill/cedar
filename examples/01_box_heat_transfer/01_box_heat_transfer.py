import cedar

# Create Mesh for Heat Transfer Model
mesh = cedar.Mesh3D("examples/01_box_heat_transfer/box_heat_transfer.vtkhdf")

# from pathlib import Path

# # Get the absolute path to the directory containing the current script
# script_dir = Path(__file__).resolve().parent

# # Define the relative path to your resource file from that directory
# relative_path_to_data = Path("data") / "input.txt"

# # Combine them to get the absolute path to the data file
# file_path = script_dir / relative_path_to_data

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
problem.add_model(ht)

problem.solve()
