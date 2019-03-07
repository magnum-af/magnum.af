import sys
import os
import arrayfire as af
import numpy as np
from magnum_af import *
import time

# expected: 
#argv[1] filepath
#argv[2] GPU number
#argv[3] nx_disk
#argv[4] a [nm]
#argv[5] Ms*mu_0
#argv[6] angle of m w.r.t. x-plane

print ("The arguments are: " , str(sys.argv))
print ("len(sys.argv)= " , len(sys.argv))
filepath = sys.argv[1]
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))

a = 1e-9*float(sys.argv[4]) if len(sys.argv) > 4 else 1.7e-6 # distance in nm between disks, argv[4] expected as a in [nm]
a_factor = a*1e6
x = 3*1e-6 + 2*a
y = 3*1e-6 + 2*a
z = 80e-9

print("a", a)
print("a_factor", a_factor)
print ("x,y,z:", x, y, z)


nx_disk = int(sys.argv[3]) if len(sys.argv) > 3 else 250
ny_disk = nx_disk
nz_disk = 1

angle_degree = float(sys.argv[6]) if len(sys.argv) > 6 else 0
angle_rad = angle_degree/360. * 2 * np.pi

disk, n_cells = Util.disk(nx_disk, ny_disk, nz_disk, axis=[np.cos(angle_rad), np.sin(angle_rad), 0])
print (n_cells)
#print (disk)

nx = int((3 + 2*a_factor) * nx_disk)
ny = int((3 + 2*a_factor) * ny_disk)
nz = nz_disk 
print (nx, ny, nz)

mesh=Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
ms_in = float(sys.argv[5])*1e-3/Constants.mu0 if len(sys.argv) > 5 else 1.75/Constants.mu0
material = Material(ms=ms_in)#TODO, confirm removal, A = 1.5e-11
m = af.constant(0., nx, ny, nz, 3, dtype=af.Dtype.f64)
m[0:nx_disk, 0:ny_disk, :, :] = disk
m[0:nx_disk, ny_disk+1+int(ny_disk*a_factor):ny_disk+1+int(ny_disk*a_factor)+ny_disk, :, :] = disk
m[0:nx_disk, -ny_disk-1:-1, :, :] = disk

m[nx_disk+1+int(nx_disk*a_factor):nx_disk+1+int(nx_disk*a_factor)+nx_disk, 0:ny_disk, :, :] = disk
#m[nx_disk+1+int(nx_disk*a_factor):nx_disk+1+int(nx_disk*a_factor)+nx_disk, ny_disk+1+int(ny_disk*a_factor):ny_disk+1+int(ny_disk*a_factor)+ny_disk, :, :] = disk
m[nx_disk+1+int(nx_disk*a_factor):nx_disk+1+int(nx_disk*a_factor)+nx_disk, -ny_disk-1:-1, :, :] = disk

m[-nx_disk-1:-1, 0:ny_disk, :, :] = disk
m[-nx_disk-1:-1, ny_disk+1+int(ny_disk*a_factor):ny_disk+1+int(ny_disk*a_factor)+ny_disk, :, :] = disk
m[-nx_disk-1:-1, -ny_disk-1:-1, :, :] = disk

state = State(mesh, material, m)
state.py_vti_writer_micro(filepath+"init")
demag = DemagField(mesh, material)
llg = LLGIntegrator(demag)
demagfield = llg.get_fheff(state)
#print (demagfield[nx/2, ny/2,:,:])
print (demagfield[nx/2, ny/2,:,0].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,1].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,2].scalar()*Constants.mu0)
stream = open(filepath+"demag.dat", "w")
stream.write("%d, %d, %e, %e, %e, %e, %e, %e" %(nx, nx_disk, a*1e9, demagfield[nx/2, ny/2,:,0].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,1].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,2].scalar()*Constants.mu0, state.material.ms, angle_degree))
stream.close()

dirty_workaround = State(mesh, material, demagfield)
dirty_workaround.py_vti_writer_micro(filepath+"demagfield")
