import sys
import os
import arrayfire as af
import numpy as np
from magnumaf import *
import time

start = time.time()
timer = time.time()
print ("The arguments are: " , str(sys.argv))
print ("len(sys.argv)= " , len(sys.argv))
filepath = sys.argv[1]
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))

x = 1e-6
y = 1e-6
z = 90e-9

print ("x, y, z:", x, y, z)


nx = 250
ny = 250
nz = 90

disk, n_cells = Util.disk(nx, ny, nz, axis=[1, 0, 0])
print (n_cells)
#print (disk)

mesh=Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
ms_in = float(sys.argv[5])*1e-3/Constants.mu0 if len(sys.argv) > 5 else 1.75/Constants.mu0
material = Material(ms=ms_in)
m = disk
m[:, :, 80:, :] = 0.
for i in range(nz):
    print(i, m[nx/2, ny/2, i, 0].scalar())

state = State(mesh, material, m)
state.write_vti(filepath+"init")
print("setup in ", time.time() - timer, "[s]")
timer = time.time()
demag = DemagField(mesh, material, verbose = True, caching = False, nthreads = 0)
print("demagtensor in ", time.time() - timer, "[s]")
timer = time.time()
llg = LLGIntegrator([demag])
demagfield = llg.h(state)
#print (demagfield[nx/2, ny/2, :, :])
print (demagfield[nx/2, ny/2, 0, 0].scalar()*Constants.mu0, demagfield[nx/2, ny/2, 0, 1].scalar()*Constants.mu0, demagfield[nx/2, ny/2, 0, 2].scalar()*Constants.mu0)

stream = open(filepath+"demag.dat", "w")
for i in range(nz):
    stream.write("%e, %e, %e\n" %(demagfield[nx/2, ny/2, i, 0].scalar()*Constants.mu0, demagfield[nx/2, ny/2, i, 1].scalar()*Constants.mu0, demagfield[nx/2, ny/2, i, 2].scalar()*Constants.mu0))
stream.close()

dirty_workaround = State(mesh, material, demagfield)
dirty_workaround.write_vti(filepath+"demagfield")
print("demagfield and output in ", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]")
