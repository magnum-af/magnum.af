#!/usr/bin/python3
# 2D DMI_D2d example from: David Cortés-Ortuño et al. Proposal for a micromagnetic standard problem for materials with Dzyaloshinskii–Moriya interaction. New J. Phys. 20, 113015 (2018).
# https://doi.org/10.1088/1367-2630/aaea1c
# https://github.com/joommf/oommf-extension-dmi-d2d

import arrayfire as af
import magnumaf as maf

args = maf.parse()  # perform CLI argument parsing (use -h flag for details)

# physical dimensions in [m]
x, y, z = 100e-9, 100e-9, 2.e-9
# discretization
nx, ny, nz = 50, 50, 1
# material parameters:
A = 13e-12  # [J/m]
D = 3.0e-3  # [J/m2]
Ms = 0.86e6  # [A/m]
Ku = 0.4e6  # [J/m3]
Ku_axis = [0.0, 0.0, 1.0]

# initial magnetization
m0 = maf.Util.disk(nx, ny, nz, axis=[0.0, 0.0, 1.0])
m0[nx/2-5: nx/2+5, ny/2-5: ny/2+5, :, 2] = -1.0

# creating magnum.af objects
mesh = maf.Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)
state = maf.State(mesh, Ms=8e5, m=m0)
state = maf.State(mesh, Ms, m=m0)
state.write_vti(args.outdir + "m_init.vti")

# define interactions
exc = maf.ExchangeField(A=1.3e-11)
ani = maf.UniaxialAnisotropyField(Ku, Ku_axis)
dmi = maf.DMI_D2d_Field(D)
terms = [exc, ani, dmi]

# setup integrator
llg = maf.LLGIntegrator(alpha=1.0, terms=terms)
outfile = open(args.outdir + "m.dat", "w")

# relaxing
while state.t < 10e-9:
    llg.step(state)
    print(state)
    print(state, file=outfile)
    if state.steps % 100 == 0:
        state.write_vti(args.outdir + "m_step" + str(state.steps) + ".vti")
state.write_vti(args.outdir + "m_final.vti")

outfile.close()
