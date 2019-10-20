# Example demonstrating the MuMAG Standard Problem 4
# Run with 'magnum.af sp4.py' or 'python3 sp4.py $PWD'

import arrayfire as af
from magnumaf import *
import sys
import time

af.info()
start = time.time()

# Physical dimensions in [m]
x, y, z =  5.e-7, 1.25e-7, 3.e-9
# Discretization
nx, ny, nz = 100, 25, 1

# Creating objects
mesh = Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)

# Initial magnetization configuration
m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
m0[1:-1, :, :, 0] = af.constant(1.0, nx-2 , ny, nz, 1, dtype=af.Dtype.f32)
m0[0, :, :, 1]    = 1.
m0[-1, :, :, 1]   = 1.

state = State(mesh, Ms = 8e5, m = m0)
demag = DemagField(mesh, verbose = True, caching = True, nthreads = 6)
exch = ExchangeField(A = 1.3e-11)
llg = LLGIntegrator(alpha = 1, terms = [demag, exch])

# Relaxing
print("relaxing 1ns")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
while state.t < 1e-9:
    llg.step(state)
    mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
    stream.write("%e, %e, %e, %e\n" %(state.t, mean[0, 0, 0, 0].scalar(), mean[0, 0, 0, 1].scalar(), mean[0, 0, 0, 2].scalar()))
print("relaxed in", time.time() - timer, "[s]")

# Preparing switch by resetting alpha and adding Zeeman field
llg.alpha=0.02
zeeswitch = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
zeeswitch[:, :, :, 0] = -24.6e-3/Constants.mu0
zeeswitch[:, :, :, 1] = +4.3e-3/Constants.mu0
zeeswitch[:, :, :, 2] = 0.0
zee = ExternalField(zeeswitch)
llg.add_terms(zee)

# Switching
print("switching 1ns")
timer = time.time()
while state.t < 2e-9:
    llg.step(state)
    mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
    stream.write("%e, %e, %e, %e\n" %(state.t, mean[0, 0, 0, 0].scalar(), mean[0, 0, 0, 1].scalar(), mean[0, 0, 0, 2].scalar()))
stream.close()
print("switched in", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]")
