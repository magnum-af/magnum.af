# Paul Heistracher <paul.thomas.heistracher@univie.ac.at>
import sys
import os
import arrayfire as af
import numpy as np
from magnum_af import *
import time
from math import fabs

# Setting filepath
print ("The arguments are: " , str(sys.argv))
filepath = sys.argv[1]
if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))
af.info()

# Physical dimensions in [m]
x = 100e-9
y = 200e-9
z = 35e-9
# Discretization 
nx = 100
ny = 100
nz = 8

## Creating mesh
mesh=Mesh(nx, ny, nz, x/nx, y/ny, z/nz)

# Setting material parameters
material=Material()
material.ms=1.58/Constants.mu0 # Saturation magnetization
material.A=15e-12 # Exchange constant
material.Ku1=0 # Anisotropy constant
material.Ku1_axis=[1, 0, 0] # TODO
material.alpha = 0.02

# Create state object with timing
start = time.time()
disk, n_cells  = Util.disk(nx, ny, nz, axis=[1,0,0])

fixed_is_elliptical=True # switch between rectangular or elliptic fixed film, set to False if needed

if fixed_is_elliptical == True:
  fixed_layer, n_cells_2  = Util.disk(nx, ny, 1, axis=[0,0,1])
  disk[:, :, 0 , :] = fixed_layer 
else:
  disk[:, :, 0 , :] = af.constant(0., nx, ny, 1, 3, dtype=af.Dtype.f64)
  disk[:, :, 0 , 2] = af.constant(1., nx, ny, 1, 1, dtype=af.Dtype.f64)


#print ("Info: n_cells = ",n_cells, " n_boolean = ", n_boolean)
state = State(mesh, material, disk)
#state.m[:, :, 0 , :] = fixed_layer # TODO does this work in wrapping? maybe this is not passed to cpp but gets stuck in wrapping
state.py_vti_writer_micro(filepath + "init_m")
print(state.meanxyz(0), state.meanxyz(1), state.meanxyz(2), np.sqrt((state.meanxyz(0))**2 +(state.meanxyz(1))**2 +(state.meanxyz(2))**2))
print ("Initialized disk configuration in ", time.time() - start, "[s]")

# Defining interaction terms
start = time.time()
fields = [
    DemagField(mesh, material, verbose = True),
    ExchangeField(mesh, material),
    #UniaxialAnisotropyField(mesh, material),
    #Zee(af.constant(0.0, nx, ny, nz, 3,dtype=af.Dtype.f64))
]
print ("Initialized interaction terms in ", time.time() - start, "[s]")
Llg = LLGIntegrator(terms = fields)

# Relaxing
print("relaxing 1ns")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
i = 0
E_diff = 1e10
E_prev = 1
while E_diff > 1e-10 and state.t < 3e-8:

  # integrate one step
  Llg.llgstep(state)

  # Fixing layer to +z direction
  temp_m = state.m
  if fixed_is_elliptical == True:
    temp_m[:, :, 0 , :] = fixed_layer 
  else:
    temp_m[:, :, 0 , :] = af.constant(0., nx, ny, 1, 3, dtype=af.Dtype.f64)
    temp_m[:, :, 0 , 2] = af.constant(1., nx, ny, 1, 1, dtype=af.Dtype.f64)
  state.m = temp_m
  
  # check energy difference every 100th step
  if i % 100 == 0:
    E_current = Llg.get_E(state)
    E_diff = fabs((E_current - E_prev)/E_current)
    print ("Step i=", i, "E_diff=", E_diff)
    E_prev = E_current

  #writing output
  stream.write("%e, %e, %e, %e\n" %(state.t, state.meanxyz(0), state.meanxyz(1), state.meanxyz(2)))
  if i % 1000 == 0:
    state.py_vti_writer_micro(filepath + "m_step_" + str(i))
  i = i +1
print("relaxed in", time.time() - timer, "[s], state.t=", state.t)

#TODO#figure out how to implement:# state.m[:, :, 0 , :] = fixed_layer
