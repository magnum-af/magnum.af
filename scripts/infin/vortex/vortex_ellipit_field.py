# Paul Heistracher <paul.thomas.heistracher@univie.ac.at>
import sys
import os
import arrayfire as af
import numpy as np
from magnumaf import *
import time

# Setting filepath
print ("The arguments are: " , str(sys.argv))
filepath = sys.argv[1]

if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))

af.info()

# Physical dimensions in [m]
x = 1600e-9
y = 1600e-9
z = 65e-9
# Discretization
nx = 250
ny = 250
nz = 1

## Creating mesh
mesh=Mesh(nx, ny, nz, x/nx, y/ny, z/nz)

# Setting material parameters
# Second material class for stress
Ms = 1.75/Constants.mu0
A = 1.5e-11
stress_Ku1=1400
stress_Ku1_axis = [1, 0, 0]

# Create state object with timing
start = time.time()
state = State(mesh, Ms = Ms, m = Util.vortex(mesh, positive_z = True))
state.write_vti(filepath + "init_m")
#print(state.micro_Ms_field)
print(state.m_mean(0), state.m_mean(1), state.m_mean(2), np.sqrt((state.m_mean(0))**2 +(state.m_mean(1))**2 +(state.m_mean(2))**2))
print ("Initialized disk configuration in ", time.time() - start, "[s]")

# Defining interaction terms
start = time.time()
terms = [
    DemagField(mesh, verbose = True, caching = True),
    ExchangeField(A),
    UniaxialAnisotropyField(stress_Ku1),
    ExternalField(af.constant(0.0, nx, ny, nz, 3,dtype=af.Dtype.f64))
]
print ("Initialized interaction terms in ", time.time() - start, "[s]")

# Creating minimizer object
minimizer = LBFGS_Minimizer(terms, tol=1e-15, maxiter=1000)

start = time.time()
minimizer.minimize(state)
state.write_vti(filepath + "m_relaxed")
print ("Relaxed initial configuration in", time.time() - start, "[s]")

# Starting minimizer loop
stream = open(filepath+"m.dat", "w")
A = 0.05/Constants.mu0 # 50mT field
steps = 100
print ("A= ", A)
for i in range(0, steps):
    phi = 2. * np.pi * i/steps;
    terms[-1].set_homogeneous_field(A * np.cos(phi), A * np.sin(phi), 0)
    start = time.time()
    minimizer.minimize(state)
    stream.write("%d, %e, %e, %e, %e, %e, %e, %e\n" %(i, state.m_mean(0), state.m_mean(1), state.m_mean(2), A * np.cos(phi), A * np.sin(phi), 0, np.sqrt((state.m_mean(0))**2 +(state.m_mean(1))**2 +(state.m_mean(2))**2)))
    stream.flush()
    print ("step ", str(i), ", phi= ", phi, ", time [s]= ", time.time() - start)
    state.write_vti(filepath + "m_"+ str(i))
stream.close()
