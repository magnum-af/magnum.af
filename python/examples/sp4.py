#!/usr/bin/python3
# Example simulating the MuMAG Standard Problem 4

import arrayfire as af
import magnumaf as maf

args = maf.parse()  # perform CLI argument parsing (use -h flag for details)

# physical dimensions in [m]
x, y, z = 500e-9, 125e-9, 3e-9

# discretization
nx, ny, nz = 100, 25, 1

# initial magnetization
m0 = af.constant(0.0, nx, ny, nz, 3, af.Dtype.f64)
m0[1:nx-1, :, :, 0] = 1.0
m0[0, :, :, 1] = 1.0
m0[-1, :, :, 1] = 1.0

# creating magnum.af objects
mesh = maf.Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)
state = maf.State(mesh, Ms=8e5, m=m0)

# define interactions
dmg = maf.DemagField(mesh)
exc = maf.ExchangeField(A=1.3e-11)

# setup integrator
llg = maf.LLGIntegrator(alpha=1.0, terms=[dmg, exc])
outfile = open(args.outdir + "m.dat", "w")

# relaxing
while state.t < 1e-9:
    llg.step(state)
    print(state, file=outfile)

# preparing switch
H_ext = af.constant(0.0, nx, ny, nz, 3, af.Dtype.f64)
H_ext[:, :, :, 0] = -24.6e-3/maf.Constants.mu0
H_ext[:, :, :, 1] = 4.3e-3/maf.Constants.mu0
ext = maf.ExternalField(H_ext)
llg.add_terms(ext)
llg.alpha = 0.02

# switching
while state.t < 2e-9:
    llg.step(state)
    print(state, file=outfile)

outfile.close()
maf.Util.plot(outputdir=args.outdir, lines=[
              'u 1:2 w l t "<mx>"', 'u 1:3 w l t "<my>"', 'u 1:4 w l t "<mz>"'])
