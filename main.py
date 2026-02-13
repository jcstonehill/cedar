import cedar
import os

mesh = cedar.Mesh3D("test.vtkhdf")

ht = cedar.HeatTransfer(mesh)
ht.T.ic = 350
ht.materials = {
    "leftvol" : cedar.materials.ConstantMaterial(1000, 1, 1),
    "rightvol" : cedar.materials.ConstantMaterial(2000, 10, 1)
}
ht.bc = {
    "left" : cedar.FixedTemperatureBC(400),
    "right" :cedar.FixedTemperatureBC(300)
}
ht.source = cedar.HeatSource(1e4)

problem = cedar.Problem([ht])
#problem.t_end = 100
problem.solve()

os.system(f"cp output.vtkhdf /mnt/c/Users/jacob/Documents/output.vtkhdf")