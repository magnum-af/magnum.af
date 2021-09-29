#!/usr/bin/python3
# Example demonstrating the MuMAG Standard Problem 4

import arrayfire as af
import magnumaf as maf
import time

start = time.time()
args = maf.parse() # perform CLI argument parsing (use -h flag for details)

# Physical dimensions in [m]
x, y, z =  5.e-7, 1.25e-7, 3.e-9
# Discretization
nx, ny, nz = 100, 25, 1

# choose numerical precision:
dtype=af.Dtype.f64

# Initial magnetization configuration
m0 = af.constant(0.0, nx, ny, nz, 3, dtype)
m0[1:nx-1, :, :, 0] = 1.0
m0[     0, :, :, 1] = 1.0
m0[    -1, :, :, 1] = 1.0

# Creating magnum.af objects
mesh = maf.Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)
state = maf.State(mesh, Ms = 8e5, m = m0)

# define interactions
dmg = maf.DemagField(mesh, verbose = True, caching = True, nthreads = 6)
exc = maf.ExchangeField(A = 1.3e-11)
fieldterms = [dmg, exc]

# setup integrator
llg = maf.LLGIntegrator(alpha = 1.0, terms = fieldterms)

# Relaxing
dt_in_sec = 1e-9
print("relaxing  for", 1e9 * dt_in_sec, "[ns] ...")
outfile = open(args.outdir + "m.dat", "w", buffering = 1)
timer = time.time()

while state.t < dt_in_sec:
    llg.step(state)
    mx, my, mz = state.mean_m()
    outfile.write("%e, %e, %e, %e\n" %(state.t, mx, my, mz))

print("relaxing  for", 1e9 * dt_in_sec, "[ns] took", time.time() - timer, "[s]")

# Preparing switch by resetting alpha and adding Zeeman field
llg.alpha=0.02
H_ext_in_Apm = [-24.6e-3/maf.Constants.mu0, 4.3e-3/maf.Constants.mu0, 0.0]

H_ext = af.constant(0.0, nx, ny, nz, 3, dtype)
H_ext[:, :, :, 0] = H_ext_in_Apm[0]
H_ext[:, :, :, 1] = H_ext_in_Apm[1]
H_ext[:, :, :, 2] = H_ext_in_Apm[2]
ext = maf.ExternalField(H_ext)
llg.add_terms(ext)

# Switching
print("switching for", 1e9 * dt_in_sec, "[ns] ...")
timer = time.time()

while state.t < 2 * dt_in_sec:
    llg.step(state)
    mx, my, mz = state.mean_m()
    outfile.write("%e, %e, %e, %e\n" %(state.t, mx, my, mz))

outfile.close()
print("switching for", 1e9 * dt_in_sec, "[ns] took", time.time() - timer, "[s]")
print("simulating   ", 2e9 * dt_in_sec, "[ns] took",  time.time() - start, "[s]")

maf.Util.plot(outputdir = args.outdir, lines = ['u ($1 * 1e9):2 w l t "<mx>"', 'u ($1 * 1e9):3 w l t "<my>"', 'u ($1 * 1e9):4 w l t "<mz>"'])
