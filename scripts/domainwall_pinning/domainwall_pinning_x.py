import arrayfire as af
from magnumaf import *
from numpy import pi
from math import sqrt
import sys
import time

if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))
af.info()
start = time.time()

# Physical dimensions in [m]
x = 8.e-8
y = 1.e-9
z = 1.e-9

# Discretization 
nx = 100
ny = 1
nz = 1

# Material
soft_Aex        = 0.25e-11
soft_ms         = 0.25 / Constants.mu0
soft_K_uni      = 1e5

hard_Aex        = 1.0e-11
hard_ms         = 1.0 / Constants.mu0
hard_K_uni      = 1e6

# Analytical Result
def H(soft_Aex, soft_ms, soft_K_uni, hard_Aex, hard_ms, hard_K_uni):
    eps_k=soft_K_uni/hard_K_uni
    eps_A=soft_Aex/hard_Aex
    eps_ms=soft_ms/hard_ms
    hard_J=hard_ms*Constants.mu0
    return 2*hard_K_uni/hard_J * ((1-eps_k*eps_A)/(1+sqrt(eps_ms*eps_A))**2)

print ("H_analytic=", H(soft_Aex, soft_ms, soft_K_uni, hard_Aex, hard_ms, hard_K_uni)*Constants.mu0, " [T]")

# Creating objects
m = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m[:nx/2,:,:,0] = af.constant( 1.0, int(nx/2) , ny, nz, 1, dtype=af.Dtype.f64)
m[:nx/2,:,:,1] = af.constant( 0.3, int(nx/2) , ny, nz, 1, dtype=af.Dtype.f64)
m[nx/2:,:,:,0] = af.constant(-1.0, int(nx/2) , ny, nz, 1, dtype=af.Dtype.f64)
m[nx/2:,:,:,1] = af.constant( 0.3, int(nx/2) , ny, nz, 1, dtype=af.Dtype.f64)

mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
material = Material(alpha=1.0, Ku1_axis=[1., 0., 0.])
state = State(mesh, material, m)
state.normalize()
state.write_vti(sys.argv[1] + "minit")
#m = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
#m[:nx/2-5,:,:,0] = af.constant(1.0, int(nx/2)-5 , ny, nz, 1, dtype=af.Dtype.f64)
#m[nx/2-5:nx/2+5,:,:,2] = af.constant(1.0, 10 , ny, nz, 1, dtype=af.Dtype.f64)
#m[nx/2+5:,:,:,0] = af.constant(-1.0, int(nx/2)-5 , ny, nz, 1, dtype=af.Dtype.f64)
#print ("int(nx/2)=",nx/2)
#print(m)

# Setting A values as field
A_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
# soft
A_field[:nx/2,:,:,:] = af.constant(soft_Aex, int(nx/2) , ny, nz, 3, dtype=af.Dtype.f64)
# hard
A_field[nx/2:,:,:,:] = af.constant(hard_Aex, int(nx/2) , ny, nz, 3, dtype=af.Dtype.f64)

# Setting Ms values as field
Ms_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
# soft
Ms_field[:nx/2,:,:,:] = af.constant(soft_ms, int(nx/2) , ny, nz, 3, dtype=af.Dtype.f64)
# hard
Ms_field[nx/2:,:,:,:] = af.constant(hard_ms, int(nx/2) , ny, nz, 3, dtype=af.Dtype.f64)

# Setting Ku1 values as field
Ku1_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
# soft
Ku1_field[:nx/2,:,:,:] = af.constant(soft_K_uni, int(nx/2) , ny, nz, 3, dtype=af.Dtype.f64)
# hard
Ku1_field[nx/2:,:,:,:] = af.constant(hard_K_uni, int(nx/2) , ny, nz, 3, dtype=af.Dtype.f64)

#print ("material:", material.alpha, state.Ms, material.A, material.Ku1, material.Ku1_axis)
#TODO# string: setting ms=0 or not setting ms leads to segfault!: #material = Material(alpha=1.0, ms=0, Ku1_axis=[1., 0., 0.])
#also not setting ms: material = Material(alpha=1.0, Ku1_axis=[1., 0., 0.])

#print (state.m)
state.micro_A_field = A_field
state.micro_Ms_field = Ms_field
state.micro_Ku1_field = Ku1_field
#state.write_vti(sys.argv[1] + "m_init")
#print (Ms_field)
#print (state.micro_Ms_field)
#print (Ku1_field)
#print (state.micro_Ku1_field)
#print (A_field)
#print (state.micro_A_field)
fields = [
    ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)),
    ExchangeField(state.mesh, state.material),
    UniaxialAnisotropyField(state.mesh, state.material),
]
Llg = LLGIntegrator(terms=fields)

#stream = open(sys.argv[1]+"m_relax.dat", "w")
#print("Relaxing 2.2e-10 [s]")
#state.write_vti(sys.argv[1] + "m_init")
#i = 0
#while state.t < 2.2e-10: #value for 100x6x6 geometry
#  if i%200 == 0:
#    state.write_vti(sys.argv[1] + "m_relax" + str(i))
#  i = i + 1
#  Llg.step(state)
#  mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
#  stream.write("%e, %e, %e, %e\n" %(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar()))
#  #stream.write("%e, %e, %e, %e, %e\n" %(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), -state.t/50e-9/Constants.mu0))
#  stream.flush()
#state.write_vti(sys.argv[1] + "m_relax")
#stream.close()

fastenup = 10
print("Start 100 [ns] hysteresis")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
i = 0
while (state.t < 1e-7/fastenup and state.m_mean(0) < (1. - 1e-6)):
  if i%2000 == 0:
    state.write_vti(sys.argv[1] + "m_" + str(i))
  fields[0].set_homogenuous_field(state, fastenup * state.t/50e-9/Constants.mu0, 0.0, 0.0)
  Llg.step(state)
  #mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
  printzee = af.mean(af.mean(af.mean(fields[0].get_zee(), dim=0), dim=1), dim=2)
  print(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), state.t/50e-9/Constants.mu0, printzee[0,0,0,0].scalar()*Constants.mu0)
  printzee = af.mean(af.mean(af.mean(fields[0].get_zee(), dim=0), dim=1), dim=2)
  stream.write("%e, %e, %e, %e, %e, %e\n" %(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), state.t/50e-9/Constants.mu0, printzee[0,0,0,0].scalar()))
  stream.flush()
  i = i + 1
print("100 [ns] hysteresis in ", time.time() - timer, "[s]")

stream.close()
print("switched at Hext=", printzee[0,0,0,0].scalar()*Constants.mu0, " [T] with state.t=", state.t, " [s]")
print("total time =", time.time() - start, "[s]")
