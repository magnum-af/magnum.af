#!/usr/bin/python3
# Example demonstrating hysteresis looping.

import arrayfire as af
from magnumaf import *
import time
import sys

# Physical dimensions in [m]
x, y, z = 1e-9, 1e-9, 1e-9
# Discretization
nx, ny, nz = 1, 1, 1

# Creating objects
mesh = Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)

# Initial magnetization configuration
m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, :, 0] = 1.

state = State(mesh, Ms = 8e5, m = m0)
external = ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64))
llg = LLGIntegrator(alpha = 1, terms = [external])

# Using external.set_homogeneous_field()
for i in range(10):
    external.set_homogeneous_field(i, i/2, 0)
    hx = external.H_in_Apm(state)[0, 0, 0, 0].scalar()
    hy = external.H_in_Apm(state)[0, 0, 0, 1].scalar()
    hz = external.H_in_Apm(state)[0, 0, 0, 2].scalar()
    print(i, hx, hy, hz)

# More general alternative
print()
for i in range(10):
    h_current = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.f64)
    h_current[:, :, :, 0] = i; # Setting Hx component
    h_current[:, :, :, 1] = i/2; # Setting Hy component
    external = ExternalField(h_current) # this overwrites the existing object named 'external'
    hx = external.H_in_Apm(state)[0, 0, 0, 0].scalar()
    hy = external.H_in_Apm(state)[0, 0, 0, 1].scalar()
    hz = external.H_in_Apm(state)[0, 0, 0, 2].scalar()
    print(i, hx, hy, hz)
