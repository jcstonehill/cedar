import gmsh
import numpy as np

f2f = 1.905/100
pitch = 0.40894/100
coating_or = 0.25654/200
flow_or = 0.23654/200

gmsh.initialize()

def hole_locations(pitch):
    locations = []

    locations.append((0, 0))
    locations.append((pitch, 0))
    locations.append((-pitch, 0))
    locations.append((2*pitch, 0))
    locations.append((-2*pitch, 0))

    x = pitch/2
    y = (pitch**2 - x**2)**0.5
    locations.append((x, y))
    locations.append((-x, y))
    locations.append((x, -y))
    locations.append((-x, -y))

    x = pitch*1.5
    locations.append((x, y))
    locations.append((-x, y))
    locations.append((x, -y))
    locations.append((-x, -y))

    x = pitch
    y = 2*y
    locations.append((x, y))
    locations.append((-x, y))
    locations.append((x, -y))
    locations.append((-x, -y))

    locations.append((0, y))
    locations.append((0, -y))

    return locations

def add_hexagon(f2f: float, cc_or: float, cc_pitch):
    r = f2f/2
    a = 2*r/np.sqrt(3)
    d = 2*a

    p1 = gmsh.model.occ.add_point(0, 0, 0)
    p2 = gmsh.model.occ.add_point(cc_or*np.cos(np.deg2rad(60)), cc_or*np.sin(np.deg2rad(60)), 0)
    p3 = gmsh.model.occ.add_point(cc_or, 0, 0)

    p4 = gmsh.model.occ.add_point(2*cc_pitch*np.cos(np.deg2rad(30)), 0, 0)
    p5 = gmsh.model.occ.add_point(2*cc_pitch*np.cos(np.deg2rad(30))-cc_or, 0, 0)
    p6 = gmsh.model.occ.add_point(2*cc_pitch*np.cos(np.deg2rad(30)), cc_or, 0)
    p7 = gmsh.model.occ.add_point(2*cc_pitch*np.cos(np.deg2rad(30))+cc_or, 0, 0)

    p8 = gmsh.model.occ.add_point(0.5*f2f, 0, 0)

    p9 = gmsh.model.occ.add_point(a*np.cos(np.deg2rad(30)), a*np.sin(np.deg2rad(30)), 0)

    p10 = gmsh.model.occ.add_point(0.5*f2f*np.cos(np.deg2rad(60)), 0.5*f2f*np.sin(np.deg2rad(60)), 0)

    p11 = gmsh.model.occ.add_point(cc_pitch*np.cos(np.deg2rad(30)), 1.5*cc_pitch, 0)
    p12 = gmsh.model.occ.add_point(cc_pitch*np.cos(np.deg2rad(30))+(np.cos(np.deg2rad(60))*cc_or), 1.5*cc_pitch+(np.sin(np.deg2rad(60))*cc_or), 0)
    p13 = gmsh.model.occ.add_point(cc_pitch*np.cos(np.deg2rad(30)), 1.5*cc_pitch-cc_or, 0)
    p14 = gmsh.model.occ.add_point(cc_pitch*np.cos(np.deg2rad(30))-(np.cos(np.deg2rad(60))*cc_or), 1.5*cc_pitch-(np.sin(np.deg2rad(60))*cc_or), 0)
    gmsh.model.occ.synchronize()


    # p4 = gmsh.model.occ.add_point(cc_pitch, 0, 0)
    # p5 = gmsh.model.occ.add_point(cc_pitch - cc_or, 0, 0)
    # p6 = gmsh.model.occ.add_point(cc_pitch, cc_or, 0)
    # p7 = gmsh.model.occ.add_point(cc_pitch + cc_or, 0, 0)

    # p8 = gmsh.model.occ.add_point(2*cc_pitch, 0, 0)
    # p9 = gmsh.model.occ.add_point(2*cc_pitch - cc_or, 0, 0)
    # p10 = gmsh.model.occ.add_point(2*cc_pitch, cc_or, 0)
    # p11 = gmsh.model.occ.add_point(2*cc_pitch + cc_or, 0, 0)

    # p12 = gmsh.model.occ.add_point(a, 0, 0)
    # p13 = gmsh.model.occ.add_point(a-(a/2)*np.cos(np.deg2rad(60)), (a/2)*np.sin(np.deg2rad(60)), 0)

    # p14 = gmsh.model.occ.add_point(1.5*pitch, np.sqrt(0.75)*pitch, 0)
    # p15 = gmsh.model.occ.add_point(1.5*pitch+cc_or*np.cos(np.deg2rad(30)), np.sqrt(0.75)*pitch+cc_or*np.sin(np.deg2rad(30)), 0)
    # p16 = gmsh.model.occ.add_point(1.5*pitch-cc_or*np.cos(np.deg2rad(120)), np.sqrt(0.75)*pitch-cc_or*np.sin(np.deg2rad(120)), 0)
    # p17 = gmsh.model.occ.add_point(1.5*pitch-cc_or*np.cos(np.deg2rad(30)), np.sqrt(0.75)*pitch-cc_or*np.sin(np.deg2rad(30)), 0)

    # gmsh.model.occ.synchronize()

    l1 = gmsh.model.occ.add_circle_arc(p2, p1, p3)
    l2 = gmsh.model.occ.add_line(p3, p5)
    l3 = gmsh.model.occ.add_circle_arc(p5, p4, p6)
    l4 = gmsh.model.occ.add_circle_arc(p6, p4, p7)
    l5 = gmsh.model.occ.add_line(p7, p8)
    l6 = gmsh.model.occ.add_line(p8, p9)
    l7 = gmsh.model.occ.add_line(p9, p10)
    l8 = gmsh.model.occ.add_line(p10, p12)
    l9 = gmsh.model.occ.add_circle_arc(p12, p11, p13)
    l10 = gmsh.model.occ.add_circle_arc(p13, p11, p14)
    l11 = gmsh.model.occ.add_line(p14, p2)

    l12 = gmsh.model.occ.add_circle(cc_pitch*np.cos(np.deg2rad(30)), cc_pitch*np.sin(np.deg2rad(30)), 0, cc_or)
    l13 = gmsh.model.occ.add_circle(2*cc_pitch*np.cos(np.deg2rad(30)), 2*cc_pitch*np.sin(np.deg2rad(30)), 0, cc_or)
    gmsh.model.occ.synchronize()
    # l6 = gmsh.model.occ.add_circle_arc(p9, p8, p10)
    # l7 = gmsh.model.occ.add_circle_arc(p10, p8, p11)
    # l8 = gmsh.model.occ.add_line(p11, p12)
    # l9 = gmsh.model.occ.add_line(p12, p13)
    # l10 = gmsh.model.occ.add_line(p13, p15)
    # l11 = gmsh.model.occ.add_circle_arc(p15, p14, p16)
    # l12 = gmsh.model.occ.add_circle_arc(p16, p14, p17)
    # l13 = gmsh.model.occ.add_line(p17, p2)
    # gmsh.model.occ.synchronize()

    C = 1
    gmsh.model.mesh.set_transfinite_curve(l1, int(C*7))
    gmsh.model.mesh.set_transfinite_curve(l2, int(C*10))
    gmsh.model.mesh.set_transfinite_curve(l3, int(C*10))
    gmsh.model.mesh.set_transfinite_curve(l4, int(C*10))
    gmsh.model.mesh.set_transfinite_curve(l5, int(C*5))
    gmsh.model.mesh.set_transfinite_curve(l6, int(C*10))
    gmsh.model.mesh.set_transfinite_curve(l7, int(C*10))
    gmsh.model.mesh.set_transfinite_curve(l8, int(C*5))
    gmsh.model.mesh.set_transfinite_curve(l9, int(C*15))
    gmsh.model.mesh.set_transfinite_curve(l10, int(C*5))
    gmsh.model.mesh.set_transfinite_curve(l11, int(C*10))
    gmsh.model.mesh.set_transfinite_curve(l12, int(C*40))
    gmsh.model.mesh.set_transfinite_curve(l13, int(C*40))


    loop1 = gmsh.model.occ.add_curve_loop([l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11])
    hole1 = gmsh.model.occ.add_curve_loop([l12])
    hole2 = gmsh.model.occ.add_curve_loop([l13])
    
    gmsh.model.occ.synchronize()
    surf = gmsh.model.occ.add_plane_surface([loop1, hole1, hole2])
    gmsh.model.occ.synchronize()

    return surf
    # p1 = gmsh.model.occ.add_point(a/2, r, 0)
    # p2 = gmsh.model.occ.add_point(a, 0, 0)
    # p3 = gmsh.model.occ.add_point(a/2, -r, 0)
    # p4 = gmsh.model.occ.add_point(-a/2, -r, 0)
    # p5 = gmsh.model.occ.add_point(-a, 0, 0)
    # p6 = gmsh.model.occ.add_point(-a/2, r, 0)
    # gmsh.model.occ.synchronize()
    # gmsh.model.occ.add_circle_arc

    # l1 = gmsh.model.occ.add_line(p1, p2)
    # l2 = gmsh.model.occ.add_line(p2, p3)
    # l3 = gmsh.model.occ.add_line(p3, p4)
    # l4 = gmsh.model.occ.add_line(p4, p5)
    # l5 = gmsh.model.occ.add_line(p5, p6)
    # l6 = gmsh.model.occ.add_line(p6, p1)
    # gmsh.model.occ.synchronize()

    # hex_loop = gmsh.model.occ.add_curve_loop([l1, l2, l3, l4, l5, l6])
    # gmsh.model.occ.synchronize()
    # loops = [hex_loop]

    # for x, y in hole_locations:
    #     circle = gmsh.model.occ.add_circle(x, y, 0, coating_or)
    #     loop = gmsh.model.occ.add_curve_loop([circle])

    #     loops.append(loop)
    #     gmsh.model.occ.synchronize()
    
    

    # surf = gmsh.model.occ.add_plane_surface(loops)
    # gmsh.model.occ.synchronize()

    # return surf

surf = add_hexagon(f2f, coating_or, pitch)

gmsh.model.occ.extrude([(2, surf)], 0, 0, 0.89, numElements = [200])
gmsh.model.occ.synchronize()
 
gmsh.model.add_physical_group(3, [1], name = "fuel")
gmsh.model.add_physical_group(2, [3, 6], name = "side1")
gmsh.model.add_physical_group(2, [7], name = "wall_to_tietube")
gmsh.model.add_physical_group(2, [8], name = "wall_to_fe")
gmsh.model.add_physical_group(2, [9, 12], name = "side2")
gmsh.model.add_physical_group(2, [1], name = "bottom")
gmsh.model.add_physical_group(2, [15], name = "top")
gmsh.model.add_physical_group(2, [2, 4, 5, 10, 11, 13, 14], name = "fc_wall")


# Set mesh size
# gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.001)
# gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.001)

# Generate the 3D mesh
gmsh.model.mesh.generate(3)

gmsh.option.set_number("Mesh.MeshSizeFactor", 1)
#gmsh.option.set_number("Mesh.SaveAll", 1)

# Write mesh file
gmsh.write("snre_fe.msh")