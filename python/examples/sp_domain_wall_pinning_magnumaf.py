#!/usr/bin/python3
# Standard problem domain wall pinning
import arrayfire as af
import  magnumaf as maf
from math import sqrt

def Hp_analytic_in_Apm(soft_Aex, soft_Ms, soft_Ku, hard_Aex, hard_Ms, hard_Ku):
    eps_k = soft_Ku / hard_Ku
    eps_A = soft_Aex / hard_Aex
    eps_ms = soft_Ms / hard_Ms
    hard_J = hard_Ms * maf.Constants.mu0
    return 2 * hard_Ku / hard_J * ((1-eps_k*eps_A)/(1+sqrt(eps_ms*eps_A))**2)

# Discretization
nx = 80
ny = nz = 1
dx = dy = dz = 1e-9

# Material parameters
soft_Aex = 0.25e-11
soft_Ku  = 1e5
soft_Ms  = 0.25 / maf.Constants.mu0

hard_Aex = 1.0e-11
hard_Ms  = 1.0 / maf.Constants.mu0
hard_Ku  = 1e6

# Initial magnetization
m0 = af.constant(0.0, nx, ny, nz, 3, af.Dtype.f64)
m0[:nx/2, :, :, 0] =  1.0
m0[nx/2:, :, :, 0] = -1.0
m0[:, :, :, 1] =  0.3
m0 = maf.Util.normalize(m0)

# Setting A, Ms and Ku values as fields
A = af.constant(soft_Aex, nx, ny, nz, 1, af.Dtype.f64)
A[nx/2:] = hard_Aex
Ms = af.constant(soft_Ms, nx, ny, nz, 1, af.Dtype.f64)
Ms[nx/2:] = hard_Ms
Ku = af.constant(soft_Ku, nx, ny, nz, 1, af.Dtype.f64)
Ku[nx/2:] = hard_Ku

# Interactions:
mesh = maf.Mesh(nx, ny, nz, dx, dy, dz)
state = maf.State(mesh, Ms, m0)
ext = maf.ExternalField(af.constant(0.0, nx, ny, nz, 3, af.Dtype.f64))
exc = maf.ExchangeField(A, mesh)
ani = maf.UniaxialAnisotropyField(Ku, [1, 0, 0])
llg = maf.LLGIntegrator(alpha=1.0, terms = [ext, exc, ani])

# Linear increasing exteral field
class IncrField:
    def __init__(self, field_rate, startField):
        self.field_rate = field_rate # [A/m/s]
        self.startField = startField # [A/m]
    def from_time(self, time):
        return self.startField + time * self.field_rate

field = IncrField(field_rate = 2./maf.Constants.mu0/100e-9, startField = 0.0)
outfile = open(maf.parse().outdir + "m.dat", "w")

while (state.mean_mx() < (1. - 1e-3)):
    ext.set_homogeneous_field(field.from_time(state.t), 0.0, 0.0)
    llg.step(state)
    outfile.write("%e, %e, %e, %e, %e\n" %(state.t, state.mean_mx(), state.mean_my(), state.mean_mz(), field.from_time(state.t)))

outfile.close()

Hp_analytic_in_T = Hp_analytic_in_Apm(soft_Aex, soft_Ms, soft_Ku, hard_Aex, hard_Ms, hard_Ku) * maf.Constants.mu0
print("Hp_analytic=  ",  Hp_analytic_in_T, " [T]")
print("Hp_calculated=", field.from_time(state.t) * maf.Constants.mu0, " [T]")
