# Example demonstrating the MuMAG Standard Problem 4
# Run with: magnum.af sp4.py

import arrayfire as af
from magnum_af import *
import sys
import time

# argv[1] is expected to be the full path to the existing output-directory
if sys.argv[1][-1] != "/":
    sys.argv[1] = sys.argv[1] + "/"
filepath = sys.argv[1]
# argv[2] is expected to be an integer setting the GPU number
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))
af.info()

start = time.time()

# Physical dimensions in [m]
x = 5.e-7
y = 1.25e-7
z = 3.e-9
# Discretization 
nx = 100
ny = 25
nz = 1

# Creating objects
mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
m = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)

material = Material()
material.ms = 8e5
material.A = 1.3e-11
material.alpha = 1.

m[1:-1,:,:,0] = af.constant(1.0, nx-2 ,ny, nz, 1, dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0, 1    ,ny, nz, 1, dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0, 1    ,ny, nz, 1, dtype=af.Dtype.f64);

state = State(mesh, material, m)
demag = DemagField(mesh, material)
exch = ExchangeField(mesh, material)
Llg = LLGIntegrator(demag, exch)

# Relaxing
print("relaxing 1ns")
stream = open(filepath+"m.dat", "w")
timer = time.time()
while state.t < 1e-9:
  Llg.llgstep(state)
  mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
  stream.write("%e, %e, %e, %e\n" %(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar()))
print("relaxed in", time.time() - timer, "[s]")

# Resetting alpha and adding Zeeman field
state.set_alpha(0.02)
zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0] = -24.6e-3/material.mu0
zeeswitch[0,0,0,1] = +4.3e-3/material.mu0
zeeswitch[0,0,0,2] = 0.0
zeeswitch = af.tile(zeeswitch, nx, ny, nz)
zee = Zee(zeeswitch)
Llg.add_terms(zee)

# Switching
print("switching 1ns")
timer = time.time()
while state.t < 2e-9:
  Llg.llgstep(state)
  mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
  stream.write("%e, %e, %e, %e\n" %(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar()))
stream.close()
print("switched in", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]")
