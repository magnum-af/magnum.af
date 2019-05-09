import sys
import os
import arrayfire as af
import numpy as np
from magnumaf import *
import time

# expected arguments:
# description [Unit] (default) 
# argv[1] filepath
# argv[2] GPU number [0-3] (0)
# argv[3] nx_disk (100)
# argv[4] a_y [nm] (1700)
# argv[5] Ms*mu_0 [mT] (1750)
# argv[6] angle [degree] of m w.r.t. x-plane (0)
# argv[7] a_x/a_y ratio [%] (100)

start = time.time()
timer = time.time()
print ("The arguments are: " , str(sys.argv))
print ("len(sys.argv)= " , len(sys.argv))
filepath = sys.argv[1]
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))

a_y = 1e-9*float(sys.argv[4]) if len(sys.argv) > 4 else 1.7e-6 # y distance in nm between disks, argv[4] expected as a_y in [nm]
a_y_factor = a_y*1e6
a_x = float(sys.argv[7])/100*a_y if len(sys.argv) > 7 else a_y # x distance in nm between disks
a_x_factor = a_x*1e6
x = 3*1e-6 + 2*a_x
y = 3*1e-6 + 2*a_y
z = 80e-9

print("a_x", a_x, "a_y", a_y)
print("a_x_factor", a_x_factor, "a_y_factor", a_y_factor)
print ("x,y,z:", x, y, z)


nx_disk = int(sys.argv[3]) if len(sys.argv) > 3 else 100
ny_disk = nx_disk
nz_disk = 1

angle_degree = float(sys.argv[6]) if len(sys.argv) > 6 else 0
angle_rad = angle_degree/360. * 2 * np.pi

disk, n_cells = Util.disk(nx_disk, ny_disk, nz_disk, axis=[np.cos(angle_rad), np.sin(angle_rad), 0])
print (n_cells)
#print (disk)

nx = int((3 + 2*a_x_factor) * nx_disk)
ny = int((3 + 2*a_y_factor) * ny_disk)
nz = nz_disk 
print (nx, ny, nz)

mesh=Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
ms_in = float(sys.argv[5])*1e-3/Constants.mu0 if len(sys.argv) > 5 else 1.75/Constants.mu0
material = Material(ms=ms_in)
m = af.constant(0., nx, ny, nz, 3, dtype=af.Dtype.f64)
m[0:nx_disk, 0:ny_disk, :, :] = disk
m[0:nx_disk, ny_disk+1+int(ny_disk*a_y_factor):ny_disk+1+int(ny_disk*a_y_factor)+ny_disk, :, :] = disk
m[0:nx_disk, -ny_disk-1:-1, :, :] = disk

m[nx_disk+1+int(nx_disk*a_x_factor):nx_disk+1+int(nx_disk*a_x_factor)+nx_disk, 0:ny_disk, :, :] = disk
#m[nx_disk+1+int(nx_disk*a_y_factor):nx_disk+1+int(nx_disk*a_y_factor)+nx_disk, ny_disk+1+int(ny_disk*a_y_factor):ny_disk+1+int(ny_disk*a_y_factor)+ny_disk, :, :] = disk
m[nx_disk+1+int(nx_disk*a_x_factor):nx_disk+1+int(nx_disk*a_x_factor)+nx_disk, -ny_disk-1:-1, :, :] = disk

m[-nx_disk-1:-1, 0:ny_disk, :, :] = disk
m[-nx_disk-1:-1, ny_disk+1+int(ny_disk*a_y_factor):ny_disk+1+int(ny_disk*a_y_factor)+ny_disk, :, :] = disk
m[-nx_disk-1:-1, -ny_disk-1:-1, :, :] = disk

state = State(mesh, material, m)
state.py_vti_writer_micro(filepath+"init")
print("setup in ", time.time() - timer, "[s]")
timer = time.time()
demag = DemagField(mesh, material, verbose = True, caching = True)
print("demagtensor in ", time.time() - timer, "[s]")
timer = time.time()
llg = LLGIntegrator([demag])
demagfield = llg.get_fheff(state)
#print (demagfield[nx/2, ny/2,:,:])
print (demagfield[nx/2, ny/2,0,0].scalar()*Constants.mu0, demagfield[nx/2, ny/2,0,1].scalar()*Constants.mu0, demagfield[nx/2, ny/2,0,2].scalar()*Constants.mu0)
stream = open(filepath+"demag.dat", "w")
stream.write("%d, %d, %e, %e, %e, %e, %e, %e, %e, %e" %(nx, nx_disk, a_y*1e9, demagfield[nx/2, ny/2,0,0].scalar()*Constants.mu0, demagfield[nx/2, ny/2,0,1].scalar()*Constants.mu0, demagfield[nx/2, ny/2,0,2].scalar()*Constants.mu0, state.material.ms, angle_degree, a_x*1e9, a_x/a_y))
stream.close()

dirty_workaround = State(mesh, material, demagfield)
dirty_workaround.py_vti_writer_micro(filepath+"demagfield")
print("demagfield and output in ", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]")
