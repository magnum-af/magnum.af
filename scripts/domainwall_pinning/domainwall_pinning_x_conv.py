import arrayfire as af
from magnumaf import *
from numpy import pi
from math import sqrt
import sys
import time

print(0.25/Constants.mu0)
print(1/Constants.mu0)
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

H_analytic = H(soft_Aex, soft_ms, soft_K_uni, hard_Aex, hard_ms, hard_K_uni)*Constants.mu0 # in [T]
# Creating objects
m = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
m[:nx/2, :, :, 0] =  1.0
m[:nx/2, :, :, 1] =  0.3
m[nx/2:, :, :, 0] = -1.0
m[nx/2:, :, :, 1] =  0.3

mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)

# Setting A values as field
A_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
# soft
A_field[:nx/2, :, :, :] = soft_Aex
# hard
A_field[nx/2:, :, :, :] = hard_Aex

# Setting Ms values as field
Ms_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
# soft
Ms_field[:nx/2, :, :, :] = soft_ms
# hard
Ms_field[nx/2:, :, :, :] = hard_ms

# Setting Ku1 values as field
Ku1_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
# soft
Ku1_field[:nx/2, :, :, :] = soft_K_uni
# hard
Ku1_field[nx/2:, :, :, :] = hard_K_uni

state = State(mesh, Ms_field, m)
state.normalize()
state.write_vti(sys.argv[1] + "minit")
fields = [
    ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)),
    ExchangeField(A_field),
    UniaxialAnisotropyField(Ku1_field, Ku1_axis=[1., 0., 0.]),
    #DemagField(mesh),
]
Llg = LLGIntegrator(alpha=1.0, terms=fields)

maxField = 2./Constants.mu0 # 2 [T]
simtime = 100e-9 # [s]

print("Start 100 [ns] hysteresis")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
i = 0
nvti = 0
while (state.t < simtime and state.m_mean(0) < (1. - 1e-6)):
    if i%300 == 0:
        state.write_vti(sys.argv[1] + "m_" + str(nvti))
        nvti = nvti + 1
    fields[0].set_homogeneous_field(state.t / simtime * maxField, 0.0, 0.0)
    Llg.step(state)
    printzee = af.mean(af.mean(af.mean(fields[0].h(state), dim=0), dim=1), dim=2)
    print(printzee[0, 0, 0, 0].scalar()*Constants.mu0 / H_analytic, state.t/simtime, state.m_mean(0), state.t /simtime * maxField, printzee[0, 0, 0, 0].scalar()*Constants.mu0)
    stream.write("%e, %e, %e, %e, %e, %e\n" %(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), state.t/50e-9/Constants.mu0, printzee[0, 0, 0, 0].scalar()))
    stream.flush()
    i = i + 1

stream.close()
state.write_vti(sys.argv[1] + "m_switched")
print("simulated", state.t, "[s] (out of", simtime, ") in", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]\n")
print("H_analytic=", H_analytic, " [T]")
print("switched at Hext=", printzee[0, 0, 0, 0].scalar()*Constants.mu0, " [T] with state.t=", state.t, " [s]")
