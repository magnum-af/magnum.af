import arrayfire as af
from magnum_af import *
from numpy import pi
import sys
import time

af.info()
start = time.time()

# Physical dimensions in [m]
x = 5.e-7
y = 5.e-9
z = 5.e-9
#y = 1.25e-7

# Discretization 
nx = 100
ny = 1
nz = 1

# Creating objects
m = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
print ("int(nx/2)=",nx/2)
m[:nx/2,:,:,0] = af.constant(1.0, int(nx/2) , ny, nz, 1, dtype=af.Dtype.f64)
m[nx/2:,:,:,0] = af.constant(-1.0, int(nx/2) , ny, nz, 1, dtype=af.Dtype.f64)
#print (m)

mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)

material = Material(alpha=1.0, ms=0.25/Constants.mu0, A=0.25e-11, Ku1=1e5, Ku1_axis=[0., 0., 1.])
print ("material:", material.alpha, material.ms, material.A, material.Ku1, material.Ku1_axis)

state = State(mesh, material, m)
demag = DemagField(mesh, material)
exch = ExchangeField(mesh, material)
aniso_field = UniaxialAnisotropyField(mesh, material)
zee_field = Zee(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64))
Llg = LLGIntegrator(demag, exch, aniso_field, zee_field)

print("Start 100 [ns] hysteresis")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
while state.t < 1e-7:
  zee_field.set_xyz(state, 0.0, 0.0, -state.t/50e-9/Constants.mu0)
  Llg.llgstep(state)
  mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
  stream.write("%e, %e, %e, %e\n" %(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar()))
print("100 [ns] hysteresis in ", time.time() - timer, "[s]")

stream.close()
print("total time =", time.time() - start, "[s]")
