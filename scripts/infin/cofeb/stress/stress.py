# Paul Heistracher <paul.thomas.heistracher@univie.ac.at>
import sys
import os
import arrayfire as af
import numpy as np
from magnumaf import *
import time

# Setting filepath
print ("The arguments are: " , str(sys.argv))
if sys.argv[1][-1] != "/":
    sys.argv[1] = sys.argv[1] + "/"
    print ("Info: Added '/' to sys.arvg[1]: ", sys.argv[1])

filepath = sys.argv[1]
#os.makedirs(filepath)# to overwrite, add: , exist_ok=True

# Initializing disk with magnetization in x, y or z
# xyz=0 initializes magnetization in x, xyz=1 in y, xyz=2 in z direction, default is 2 == z
def disk(n0, n1, n2, xyz = 2):
    n_cells=0
    m = np.zeros((n0, n1, n2, 3));
    for ix in range (0, n0):
        for iy in range(0, n1):
            for iz in range(0, n2):
                a= n0/2
                b= n1/2
                rx=ix-n0/2.
                ry=iy-n1/2.
                r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2);
                if(r<1):
                    m[ix, iy, iz, xyz]=1
                    n_cells = n_cells +1
    return af.from_ndarray(m), n_cells

# Initializing boolean array where only values with 1 taken into account in the calculation of the mean magnetization
def boolean_disk(n0, n1, n2, r_inner = 0.9):
    n_cells=0
    boolean = np.zeros((n0, n1, n2), dtype = bool);
    for ix in range (0, n0):
        for iy in range(0, n1):
            for iz in range(0, n2):
                a= n0/2
                b= n1/2
                rx=ix-n0/2.
                ry=iy-n1/2.
                r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2);
                if(r < r_inner):# NOTE: (keep in mind that in general 'r' is not the radius of a circle and for e.g. r2=2*r1, A2 != 4*A1)
                    boolean[ix, iy, iz]=1
                    n_cells = n_cells +1
    return af.from_ndarray(boolean), n_cells

# Setting GPU number (0-3 aviable on GTO)
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))

#af.set_backend("cpu")# TODO this segfaults
#af.set_backend("opencl")#TODO this segfaults on GTO
af.info()

# Physical dimensions in [m]
x = 800e-9
y = 800e-9
z = 1.3e-3/1.056e6
# Discretization
nx = 250
ny = 250
nz = 1

## Creating mesh
mesh=Mesh(nx, ny, nz, x/nx, y/ny, z/nz)

# Setting material parameters
material=Material()
state.Ms=1.58/Constants.mu0 # Saturation magnetization
material.A=15e-12 # Exchange constant
material.Ku1=1.3e-3/z # Anisotropy constant

# Second material class for stress
param_stress=Material()
param_stress.ms=1.58/Constants.mu0
param_stress.A=15e-12
param_stress.Ku1=1400  #TODO guessed worst case value fom Toni, elaborate
param_stress.Ku1_axis=[1, 0, 0] # Setting axis in x-direction
#print ("Check: Ku1 axis =", param_stress.print_Ku1_axis())

# Create state object with timing
start = time.time()
disk1, n_cells  = disk(nx, ny, nz)
boolean, n_boolean  = boolean_disk(nx, ny, nz, 0.9) # TODO: add respective value here

state = State(mesh, material, disk1, boolean)# NOTE update: optional argument 'boolean' allows for specified mean value evaluations
state.write_vti(filepath + "init_m")
state.write_vti_boolean(filepath + "boolean")
print(state.mean_m(0), state.mean_m(1), state.mean_m(2), np.sqrt((state.mean_m(0))**2 +(state.mean_m(1))**2 +(state.mean_m(2))**2))
print ("Info: n_cells = ", n_cells, " n_boolean = ", n_boolean)
print ("Initialized disk configuration in ", time.time() - start, "[s]")

# Defining interaction terms
start = time.time()
demag = DemagField(mesh, material)
exch=ExchangeField(mesh, material)
aniso_z = UniaxialAnisotropyField(mesh, material)
aniso_stress = UniaxialAnisotropyField(mesh, param_stress)
zee = ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64))
print ("Initialized interaction terms in ", time.time() - start, "[s]")

# Creating minimizer object
minimizer = LBFGS_Minimizer(terms=[demag, exch, aniso_z, aniso_stress, zee], tol=1e-15, maxiter=1000)

# Starting minimizer loop
stream = open(filepath+"m.dat", "w")
A = 0.05/Constants.mu0
steps = 100
print ("A= ", A)
for i in range(0, steps):
    phi = 2. * np.pi * i/steps;
    zee.set_homogeneous_field(state, A * np.cos(phi), A * np.sin(phi), 0)
    start = time.time()
    minimizer.minimize(state)
    stream.write("%d, %e, %e, %e, %e, %e, %e, %e\n" %(i, state.mean_m(0), state.mean_m(1), state.mean_m(2), A * np.cos(phi), A * np.sin(phi), 0, np.sqrt((state.mean_m(0))**2 +(state.mean_m(1))**2 +(state.mean_m(2))**2)))
    stream.flush()
    print ("step ", str(i), ", phi= ", phi, ", time [s]= ", time.time() - start)
    state.write_vti(filepath + "m_"+ str(i))
stream.close()
