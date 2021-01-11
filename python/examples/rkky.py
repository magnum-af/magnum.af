#!/usr/bin/python3
# RKKY example from https://mumax.github.io/examples.html

import numpy as np
import arrayfire as af
from magnumaf import *
import sys
import time

args = parse()
af.info()

# Discretization
nx, ny, nz = 10, 10, 2
dx, dy, dz = 1e-9, 1e-9, 1e-9

# Material parameters
Ms = 1e6
A = 10e-12
RKKY = -1e-3 * dz # TODO move to source

# Initial magnetization configuration
m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, :, 0] = 1.

# Creating objects
mesh = Mesh(nx, ny, nz, dx, dy, dz)
state = State(mesh, Ms, m = m0)

demag = DemagField(mesh, verbose = True, caching = False, nthreads = 6)

rkky_values = af.constant(RKKY/2., nx, ny, nz, 1, dtype=af.Dtype.f64)
exch_values = af.constant(A, nx, ny, nz, 1, dtype=af.Dtype.f64)
rkky_indices = af.constant(0, nx, ny, nz, 1, dtype=af.Dtype.u32)
rkkyexch = RKKYExchangeField(rkky_values, exch_values, mesh, rkky_indices)

llg = LLGIntegrator(alpha = 1, terms = [demag, rkkyexch])

E = []
print("Start rotating")
stream = open(args.outdir + "m.dat", "w")
timer = time.time()
for i in range(0, 360):
    mix = np.cos(i * np.pi/180.);
    miy = np.sin(i * np.pi/180.);
    m = state.m
    m[:, :, 1, 0] = mix
    m[:, :, 1, 1] = miy
    state.m = m
    E.append( llg.Eeff_in_J(state) )
    print("angle=[°]", i, " E=", E[-1])
    mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
    stream.write("%d, %e, %e, %e, %e\n" %(i, E[-1], mean[0, 0, 0, 0].scalar(), mean[0, 0, 0, 1].scalar(), mean[0, 0, 0, 2].scalar()))
print("fullrotation in ", time.time() - timer, "[s]")
stream.close()

print("E diff =", max(E) - min(E), "[J], should be around 2.1e-19 (from mumax3 plot)")

# plotting data with gnuplot
from os import system
system('gnuplot -e "\
    set terminal pdf;\
    set output \'' + args.outdir + 'm.pdf\';\
    set xlabel \'angle [°]\';\
    set ylabel \'E\';\
    p \'' + args.outdir + '/m.dat\' u 1:2 w l t \'E\';\
    set ylabel \'<m>\';\
    p \'' + args.outdir + '/m.dat\' u 1:3 w l t \'<m_x>\',\
    \'\' u 1:4 w l t \'<m_y>\',\
    \'\' u 1:5 w l t \'<m_z>\';\
"')

# show pdf with evince
system('evince ' + args.outdir +'m.pdf')
