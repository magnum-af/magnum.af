# Paul Heistracher <paul.thomas.heistracher@univie.ac.at>
import sys
import os
import arrayfire as af
import numpy as np
from magnum_af import *
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
            a= n0/2
            b= n1/2
            rx=ix-n0/2.
            ry=iy-n1/2.
            r = pow(rx,2)/pow(a,2)+pow(ry,2)/pow(b,2);
            if(r<1):
                m[ix,iy,:,xyz]=1
                n_cells = n_cells +1
    return af.from_ndarray(m), n_cells

# Setting GPU number (0-3 aviable on GTO)
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))
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
mesh=pyMesh(nx, ny, nz, x/nx, y/ny, z/nz)

# Setting material parameters
param=pyParam()
param.ms(1.58/param.print_mu0()) # Saturation magnetization
param.A(15e-12) # Exchange constant
param.Ku1(1.3e-3/z) # Anisotropy constant

# Second param class for stress
param_stress=pyParam()
param_stress.ms(1.58/param.print_mu0())
param_stress.A(15e-12)
param_stress.Ku1(1400) #TODO guessed worst case value fom Toni, elaborate
param_stress.Ku1_axis(1, 0, 0) # Setting axis in x-direction
print ("Check: Ku1 axis =", param_stress.print_Ku1_axis())

# Create state object with timing
start = time.time()
disk1, n_cells  = disk(nx, ny, nz)
x_corrector = (nx * ny * nz)/n_cells
y_corrector = (nx * ny * nz)/n_cells
z_corrector = (nx * ny * nz)/n_cells
state = pyState(mesh, param, disk1)
print ("n_cells",n_cells)
print ("Initialize disk configuration [s]= ", time.time() - start)

# Defining interaction terms
start = time.time()
demag = pyDemagSolver(mesh, param)
exch=pyExchSolver(mesh, param)
aniso_z = pyMicroAniso(mesh, param)
aniso_stress = pyMicroAniso(mesh, param_stress)
zee = pyZee(af.constant(0.0, nx, ny, nz, 3,dtype=af.Dtype.f64))
print ("Init terms [s]= ", time.time() - start)

# Creating minimizer object
minimizer = pyLbfgsMinimizer(demag, exch, aniso_z, aniso_stress, zee)

# Starting minimizer loop
stream = open(filepath+"m.dat", "w")
A = 0.05/param.print_mu0()
steps = 100
print ("A= ", A)
for i in range(0, steps):
    phi = 2. * np.pi * i/steps;
    zee.set_xyz(state, A * np.cos(phi), A * np.sin(phi), 0)
    start = time.time()
    minimizer.pyMinimize(state)
    stream.write("%d, %e, %e, %e, %e\n" %(i, state.meanxyz(0)*x_corrector, state.meanxyz(1)*y_corrector, state.meanxyz(2)*z_corrector, np.sqrt((state.meanxyz(0)*x_corrector)**2 +(state.meanxyz(1)*y_corrector)**2 +(state.meanxyz(2)*z_corrector)**2)))
    stream.flush()
    print ("step ", str(i), ", phi= ", phi, ", time [s]= ", time.time() - start)
    state.py_vti_writer_micro(filepath + "m_"+ str(i))
stream.close()
