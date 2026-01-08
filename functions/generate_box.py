import gmsh

# Initialize gmsh
gmsh.initialize()
gmsh.model.add("box_mesh")

# Box dimensions
length = 1.0
width = 0.5
height = 2.

# Create the box
box = gmsh.model.occ.addBox(0, 0, 0, length, width, height)

# Synchronize CAD kernel
gmsh.model.occ.synchronize()

# Get all surfaces (faces) of the box
surfaces = gmsh.model.getEntities(dim=2)

# Create physical groups (with descriptive names)
gmsh.model.addPhysicalGroup(2, [3], name="bottom")
gmsh.model.addPhysicalGroup(2, [4], name="top")
gmsh.model.addPhysicalGroup(2, [1], name="left")
gmsh.model.addPhysicalGroup(2, [2], name="right")
gmsh.model.addPhysicalGroup(2, [6], name="front")
gmsh.model.addPhysicalGroup(2, [5], name="back")
# gmsh.model.occ.synchronize()

# Add the volume as a physical group too
volumes = gmsh.model.getEntities(dim=3)
gmsh.model.addPhysicalGroup(3, [volumes[0][1]], name="volume")

# Set mesh size
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1)

# Generate the 3D mesh
gmsh.model.mesh.generate(3)

gmsh.option.set_number("Mesh.MeshSizeFactor", 1)
#gmsh.option.set_number("Mesh.MshFileVersion", 2.2)
#gmsh.option.set_number("Mesh.SaveAll", 1)

# Write mesh file
gmsh.write("box_test.msh")

# # Optionally visualize
# if "-nopopup" not in sys.argv:
#     gmsh.fltk.run()

gmsh.finalize()