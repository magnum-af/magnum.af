#!/usr/bin/python3

import arrayfire as af
from magnumaf import *
import sys
import time

args = parse()

# Discretization
nx, ny, nz = 1, 1, 1
dxyz = 1e-9

Ku = 7e6 * dxyz**3 # atom uniaxial aniso in [J]
alpha = 0.1 # damping
T = 10 # temperature in [K]
dt = 1e-15 # time step in [s]

# Initial magnetization configuration
m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, :, 0]    = 1. # set in x-dir

# External field
ext_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
ext_field[:, :, :, 0] = 20e-3/Constants.mu0 # 20mT in x

# Creating af objects
mesh = Mesh(nx, ny, nz, dx=dxyz, dy=dxyz, dz=dxyz)
state = State(mesh, Ms = 8e5, m = m0)
demag = AtomisticDipoleDipoleField(mesh)
exch = AtomisticUniaxialAnisotropyField(K_atom = 1e3, K_atom_axis = [1, 0, 0])
ext = AtomisticExternalField(ext_field)

# Choose between stochastic Heun or SemiImplicitHeun method:
if True:
    smode = "Heun"
else:
    smode = "SmiHeun"

sllg = Stochastic_LLG(alpha = alpha, T = T, dt = dt, state = state, terms = [demag, exch, ext], smode = smode)
for i in range(10):
    print(i, state.mean_m())
    sllg.step(state)
