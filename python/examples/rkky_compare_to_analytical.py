#!/usr/bin/python3
#  RKKY-coupled layers with comparison to analytical solution

from os import system
import numpy as np
import arrayfire as af
from magnumaf import *
import sys
import time

args = parse()
af.info()

minimize = False  # True uses minimizer, False uses llg integration
hys_steps = 200  # hysteresis steps

# Discretization
nx, ny, nz = 10, 10, 2
x, y, z = 1e-8, 1e-8, 2e-9
dx, dy, dz = x/nx, y/ny, z/nz
print("dx={}, dy={}, dz={}".format(dx, dy, dz))


def m_analytical_nocap(dz, Ms, H, RKKY):
    return (dz * Ms * H)/(4 * abs(RKKY))


def m_analytical(dz, Ms, H, RKKY):
    result = m_analytical_nocap(dz, Ms, H, RKKY)
    if result > 1.:
        return 1
    elif result < -1.:
        return -1
    else:
        return result


def H_analytical(RKKY_surface, dz, Ms):
    return (4. * abs(RKKY_surface))/(dz * Ms)


# Material parameters
Js = 1  # [T]
Ms = Js/Constants.mu0  # "[T/mu0]"
A = 13e-12
RKKY_surface = -0.1e-3
RKKY = RKKY_surface * dz
analytical_max_field = H_analytical(RKKY_surface, dz, Ms)
print('Analytical max field [T]  =', analytical_max_field)
print('Analytical max field [A/m]=', analytical_max_field/Constants.mu0)
ext_max_field = 1.1 * analytical_max_field/Constants.mu0

# Creating objects
mesh = Mesh(nx, ny, nz, dx, dy, dz)

m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, 0, 2] = -1.
m0[:, :, 1, 2] = 1.
state = State(mesh, Ms, m=m0)
state.write_vti(args.outdir + "m_init")


RKKYarr = af.constant(RKKY, nx, ny, nz, 1, dtype=af.Dtype.f64)
excharr = af.constant(0, nx, ny, nz, 1, dtype=af.Dtype.f64)
rkkyexch = RKKYExchangeField(RKKYarr, excharr, mesh, rkky_indices=af.constant(
    0, nx, ny, nz, 1, dtype=af.Dtype.u32))


ext_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
ext = ExternalField(ext_field)
terms = [rkkyexch, ext]

if minimize:
    minimizer = LBFGS_Minimizer(terms, tol=1e-15, maxiter=1000)
else:
    # Note: for this example hysteresis 'ears' are smaller for hmax = 1e-9 is smaller than for 3.5e-10, interestingly
    llg = LLGIntegrator(alpha=1, terms=terms, hmax=1e-9)


def hysteresis_factor(i, steps):
    if i < steps/4.:
        return 4 * i/steps
    elif i < 3 * steps/4.:
        return (4 * -i/steps + 2.)
    else:
        return (4 * i/steps - 4.)


# running hysteresis loop
stream = open(args.outdir + "m.dat", "w", buffering=1)
stream.write("# Hext [T], mx, my, mz")
for i in range(1, hys_steps + 1):
    # if i > hys_steps/4:
    #    break
    extfield = hysteresis_factor(i, hys_steps) * ext_max_field
    ext.set_homogeneous_field(extfield, 0, 0)
    if minimize:
        minimizer.minimize(state)
    else:
        llg.relax(state, precision=1e-6, verbose=False)
    state.write_vti(args.outdir + "m_step_" + str(i))
    mx, my, mz = state.mean_m()
    m_a = m_analytical(dz, Ms, extfield * Constants.mu0, RKKY_surface)
    print(i, 'ext[T]={:2.3f}, mx={:1.3f}, my={:1.3f}, mz={:1.3f}, ma={:1.3f}'.format(
        ext.H_in_Apm(state)[0, 0, 0, 0].scalar() * Constants.mu0, mx, my, mz, m_a))
    stream.write("%e, %e, %e, %e, %e\n" %
                 (extfield * Constants.mu0, mx, my, mz, m_a))

stream.close()

# plotting data with gnuplot
system('gnuplot -e "\
    set terminal pdf;\
    set output \'' + args.outdir + 'm.pdf\';\
    set xlabel \'Hx [ns]\';\
    set ylabel \'<m>\';\
    p \'' + args.outdir + '/m.dat\' u 1:2 w l t \'<m_x>\',\
    \'\' u 1:5 w l t \'<m_{analytical}>\';\
"')

# show pdf with evince
system('evince ' + args.outdir + 'm.pdf')
