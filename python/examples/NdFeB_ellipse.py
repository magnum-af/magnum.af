#!/usr/bin/python3
# NdFeB-ellipse hysteresis

import arrayfire as af
from magnumaf import *
import numpy as np
import time
start_time = time.time()

args = parse() # parses comand-line arguments via sys.argv
af.info()
dtype = af.Dtype.f64 # setting precision

# Defining material parameters
Js = 1.750             # spontaneous polarization [J]
Ms = Js/Constants.mu0  # saturation magnetnetization [J/T/m^3]
A = 15.e-12            # exchange constant [J/m]
alpha = 0.02           # damping constant []
Ku1 = 20000            # anisotropy constant [J/m^3]
Ku1_axis = [1, 0, 0]

x, y, z =  500e-9, 500e-9, 14e-9 # Physical dimensions in [m]
nx, ny, nz = 64, 64, 2 # Number of cells per axis
dx, dy, dz = x/nx, y/ny, z/nz # discretization:

RKKY_surface = -3e-3 # [J/m^2]
RKKY = RKKY_surface * dz

# Switch between LLG and LBFGS Minimizer: True is LLG, False is Minimizer
llg_over_minimizer = True

# Defining geometry
geom = Geometry.xy_ellipse(nx, ny, nz, make_3d = True)

# Initial magnetization:
if True:
    # creating isotropic random distribution of unit spins in geometry
    m0 = geom * Magnetization.isotropic(nx, ny, nz, dtype)
else:
    # Alternative: ellipse with bottom +x, top -x
    m0 = af.constant(0.0, nx, ny, nz, 3, dtype)
    m0[:, :, 0, 0] = 1.  # setting mx in lower plane to 1
    m0[:, :, 1, 0] = -1. # setting mx in upper plane to -1
    m0 = geom * m0

# Creating magnumaf objects:
mesh = Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)
state = State(mesh, Ms, m = m0)
print("dx, dy, dz=", mesh.dx, mesh.dy, mesh.dz)
state.write_vti(args.outdir + "m0")

demag = DemagField(mesh, verbose = True, caching = True, nthreads = 8)
RKKY_arr = af.constant(RKKY, nx, ny, nz, 1, dtype)
exch_arr = af.constant(A, nx, ny, nz, 1, dtype)
rkkyexch = RKKYExchangeField(RKKY_arr, exch_arr, mesh)
aniso = UniaxialAnisotropyField(Ku1, Ku1_axis)
ext = ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype))
terms = [demag, rkkyexch, aniso, ext]

if llg_over_minimizer:
    llg = LLGIntegrator(alpha, terms, dissipation_term_only = True)
else:
    minimizer = LBFGS_Minimizer(terms)

# Defining Hysteresis parameters
amplitude = 0.2/Constants.mu0 # Amplitude in [A/m]
periods = 3
rad_per_step = 2.0 * np.pi / 1000.0
steps = int(periods * np.pi * 2.0 / rad_per_step)

vti_every = 10
stream = open(args.outdir + "m.dat", "w")

print("Starting hysteresis with", steps, "steps.")
for i in range(steps):
    # make step:
    ext.set_homogeneous_field(amplitude * np.sin(i * rad_per_step), 0, 0)
    if llg_over_minimizer is True:
        llg.relax(state, precision = 1e-4, verbose = False)
    else:
        minimizer.minimize(state)

    # write output:
    mx, my, mz = state.mean_m()
    Hx, Hy, Hz = Util.spacial_mean(ext.H_in_T(state))
    print(i, mx, my, mz, Hx, Hy, Hz)
    stream.write("%e, %e, %e, %e, %e, %e, %e\n" %(i, mx, my, mz, Hx, Hy, Hz))
    if steps % vti_every == 0:
        state.write_vti(args.outdir + "step" + str(i))

print("total time =", time.time() - start_time, "[s]")

Util.plot(outputdir = args.outdir, datafile = "m.dat", xlabel = 'Hx', ylabel = "<mx>", lines = ['u 5:2 w l'])
