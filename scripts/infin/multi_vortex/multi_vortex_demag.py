import sys
import os
import arrayfire as af
import numpy as np
from magnum_af import *
import time

print ("The arguments are: " , str(sys.argv))
filepath = sys.argv[1]
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))

#a = 1.7e-6
#a_factor = 1.7
x = 3*1e-6 + 2*1.7e-6
y = 3*1e-6 + 2*1.7e-6
z = 80e-9


nx_disk = int(sys.argv[3]) if len(sys.argv) > 3 else 250
ny_disk = nx_disk
nz_disk = 1

disk, n_cells = Util.disk(nx_disk, ny_disk, nz_disk, 0)
print (n_cells)
#print (disk)

nx = int((3 + 2*1.7) * nx_disk)
ny = int((3 + 2*1.7) * ny_disk)
nz = nz_disk 
print (nx, ny, nz)

mesh=Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
material = Material(ms=1.75/Constants.mu0, A = 1.5e-11)
m = af.constant(0., nx, ny, nz, 3, dtype=af.Dtype.f64)
m[0:nx_disk, 0:ny_disk, :, :] = disk
m[0:nx_disk, ny_disk+1+int(ny_disk*1.7):ny_disk+1+int(ny_disk*1.7)+ny_disk, :, :] = disk
m[0:nx_disk, -ny_disk-1:-1, :, :] = disk

m[nx_disk+1+int(nx_disk*1.7):nx_disk+1+int(nx_disk*1.7)+nx_disk, 0:ny_disk, :, :] = disk
m[nx_disk+1+int(nx_disk*1.7):nx_disk+1+int(nx_disk*1.7)+nx_disk, ny_disk+1+int(ny_disk*1.7):ny_disk+1+int(ny_disk*1.7)+ny_disk, :, :] = disk
m[nx_disk+1+int(nx_disk*1.7):nx_disk+1+int(nx_disk*1.7)+nx_disk, -ny_disk-1:-1, :, :] = disk

m[-nx_disk-1:-1, 0:ny_disk, :, :] = disk
m[-nx_disk-1:-1, ny_disk+1+int(ny_disk*1.7):ny_disk+1+int(ny_disk*1.7)+ny_disk, :, :] = disk
m[-nx_disk-1:-1, -ny_disk-1:-1, :, :] = disk

state = State(mesh, material, m)
state.py_vti_writer_micro(filepath+"init")
demag = DemagField(mesh, material)
llg = LLGIntegrator(demag)
demagfield = llg.get_fheff(state)
#print (demagfield[nx/2, ny/2,:,:])
print (demagfield[nx/2, ny/2,:,0].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,1].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,2].scalar()*Constants.mu0)
stream = open(filepath+"demag.dat", "w")
stream.write("%d, %d, %e, %e, %e" %(nx, nx_disk, demagfield[nx/2, ny/2,:,0].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,1].scalar()*Constants.mu0, demagfield[nx/2, ny/2,:,2].scalar()*Constants.mu0))
stream.close()
