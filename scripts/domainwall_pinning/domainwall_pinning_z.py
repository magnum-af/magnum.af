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
x = 1.e-9
y = 1.e-9
z = 8.e-8

# Discretization 
nx = 1
ny = 1
nz = 50

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
m[:,:,:nz/2,2] = af.constant( 1.0, nx, ny, int(nz/2), 1, dtype=af.Dtype.f64)
m[:,:,:nz/2,1] = af.constant( 0.3, nx, ny, int(nz/2), 1, dtype=af.Dtype.f64)
m[:,:,nz/2:,2] = af.constant(-1.0, nx, ny, int(nz/2), 1, dtype=af.Dtype.f64)
m[:,:,nz/2:,1] = af.constant( 0.3, nx, ny, int(nz/2), 1, dtype=af.Dtype.f64)

mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)

# Setting A values as field
A_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
# soft
A_field[:,:,:nz/2,:] = af.constant(soft_Aex, nx, ny, int(nz/2), 3, dtype=af.Dtype.f64)
# hard
A_field[:,:,nz/2:,:] = af.constant(hard_Aex, nx, ny, int(nz/2), 3, dtype=af.Dtype.f64)

# Setting Ms values as field
Ms_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
# soft
Ms_field[:,:,:nz/2,:] = af.constant(soft_ms, nx, ny, int(nz/2), 3, dtype=af.Dtype.f64)
# hard               
Ms_field[:,:,nz/2:,:] = af.constant(hard_ms, nx, ny, int(nz/2), 3, dtype=af.Dtype.f64)

# Setting Ku1 values as field
Ku1_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
# soft
Ku1_field[:,:,:nz/2,:] = af.constant(soft_K_uni, nx, ny, int(nz/2), 3, dtype=af.Dtype.f64)
# hard                
Ku1_field[:,:,nz/2:,:] = af.constant(hard_K_uni, nx, ny, int(nz/2), 3, dtype=af.Dtype.f64)

state = State(mesh, Ms_field, m)
state.normalize()
state.write_vti(sys.argv[1] + "minit")

fields = [
    ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)),
    ExchangeField(A_field),
    UniaxialAnisotropyField(Ku1_field),
]
Llg = LLGIntegrator(alpha=1.0, terms=fields)

fastenup = 10
print("Start [ns] hysteresis")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
i = 0

printzee = af.mean(af.mean(af.mean(fields[0].h(state), dim=0), dim=1), dim=2)
#print("m", state.m)
#print(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), fastenup * state.t/50e-9/Constants.mu0, printzee[0,0,0,0].scalar())
#while i < 10:
while (state.t < 1e-7/fastenup and state.m_mean(2) < (1. - 1e-6)):
  #print("m", state.m)
  if i%2000 == 0:
    state.write_vti(sys.argv[1] + "m_" + str(i))
  fields[0].set_homogenuous_field(0.0, 0.0, fastenup * state.t/50e-9/Constants.mu0)
  Llg.step(state)
  printzee = af.mean(af.mean(af.mean(fields[0].h(state), dim=0), dim=1), dim=2)
  print(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), fastenup * state.t/50e-9/Constants.mu0, printzee[0,0,0,2].scalar()*Constants.mu0)
  stream.write("%e, %e, %e, %e, %e, %e\n" %(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), fastenup * state.t/50e-9/Constants.mu0, printzee[0,0,0,2].scalar()))
  stream.flush()
  i = i + 1
print("hysteresis for state.t=", state.t, " [s] in ", time.time() - timer, "[s]")
stream.close()
print("switched at Hext=", printzee[0,0,0,2].scalar()*Constants.mu0, " [T] with state.t=", state.t, " [s]")
print("total time =", time.time() - start, "[s]")
