#!/usr/bin/python3
import arrayfire as af
from magnumaf import *
from math import sqrt
import sys
import time

args = parse()
af.info()
start = time.time()

# setting wire direction (x, y or z)
print(args.posargs)
wire_dir = args.posargs[0] if len(args.posargs) > 0 else "x"
print("Setting wire direction to", wire_dir)

nn = 80 # number of cells in wire_dir
# Physical dimensions in [m]
if wire_dir == "z":
    x = 1.e-9
    y = 1.e-9
    z = 8.e-8
    nx = 1
    ny = 1
    nz = nn
    wire_dir_val = 2
elif wire_dir == "y":
    x = 1.e-9
    y = 8.e-8
    z = 1.e-9
    nx = 1
    ny = nn
    nz = 1
    wire_dir_val = 1
elif wire_dir == "x":
    x = 8.e-8
    y = 1.e-9
    z = 1.e-9
    nx = nn
    ny = 1
    nz = 1
    wire_dir_val = 0
else:
    raise ValueError("Unvalid wire dir, pease choose either {x, y, z}")

print("discretization:", x, y, z, nx, ny, nz)

# Discretization

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
print ("H_analytic=", H_analytic, " [T]")

# setting initial m (maybe rotation should be considered)
m = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
if wire_dir == "z":
    m[:, :, :nz/2, 2] =  1.0
    m[:, :, :nz/2, 1] =  0.3
    m[:, :, nz/2:, 2] = -1.0
    m[:, :, nz/2:, 1] =  0.3
elif wire_dir == "y":
    m[:, :ny/2, :, 1] =  1.0
    m[:, :ny/2, :, 2] =  0.3
    m[:, ny/2:, :, 1] = -1.0
    m[:, ny/2:, :, 2] =  0.3
elif wire_dir == "x":
    m[:nx/2, :, :, 0] =  1.0
    m[:nx/2, :, :, 1] =  0.3
    m[nx/2:, :, :, 0] = -1.0
    m[nx/2:, :, :, 1] =  0.3

m = Util.normalize(m)

mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)

# Setting A, Ms and Ku1 values as fields
A_field = af.constant(0.0, nx, ny, nz, 1, dtype=af.Dtype.f64)
Ms_field = af.constant(0.0, nx, ny, nz, 1, dtype=af.Dtype.f64)
Ku1_field = af.constant(0.0, nx, ny, nz, 1, dtype=af.Dtype.f64)

if wire_dir == "z":
    A_field  [:, :, :nz/2] = soft_Aex # soft
    A_field  [:, :, nz/2:] = hard_Aex # hard
    Ms_field [:, :, :nz/2] = soft_ms # soft
    Ms_field [:, :, nz/2:] = hard_ms # hard
    Ku1_field[:, :, :nz/2] = soft_K_uni # soft
    Ku1_field[:, :, nz/2:] = hard_K_uni # hard
    Ku1_axis=[0, 0, 1]
elif wire_dir == "y":
    A_field  [:, :ny/2] = soft_Aex # soft
    A_field  [:, ny/2:] = hard_Aex # hard
    Ms_field [:, :ny/2] = soft_ms # soft
    Ms_field [:, ny/2:] = hard_ms # hard
    Ku1_field[:, :ny/2] = soft_K_uni # soft
    Ku1_field[:, ny/2:] = hard_K_uni # hard
    Ku1_axis=[0, 1, 0]
elif wire_dir == "x":
    A_field  [:nx/2] = soft_Aex # soft
    A_field  [nx/2:] = hard_Aex # hard
    Ms_field [:nx/2] = soft_ms # soft
    Ms_field [nx/2:] = hard_ms # hard
    Ku1_field[:nx/2] = soft_K_uni # soft
    Ku1_field[nx/2:] = hard_K_uni # hard
    Ku1_axis=[1, 0, 0]

state = State(mesh, Ms_field, m)
state.write_vti(args.outdir + "minit")

fields = [
    ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)),
    SparseExchangeField(A_field, mesh),
    UniaxialAnisotropyField(Ku1_field, Ku1_axis),
]

##TODO# maxField = 0. # Note: for zero field, the domainwall does not pin to the interface for current implementation i.e. sparse yielding correct Hanalytic
# TODO minimizer not working and yielding mean_m = nan

llg_before_mini = True
if llg_before_mini == True:
    llg = LLGIntegrator(alpha=1.0, terms=fields)
else:
    minimizer = LBFGS_Minimizer(fields)

class Field:
    def __init__(self, maxField = 2./Constants.mu0 , simtime = 100e-9, startField = 0/Constants.mu0):
        self.maxField = maxField # [Oe]
        self.simtime = simtime # [s]
        self.startField = startField # [Oe]
    def from_time(self, time):
        return self.startField + time / self.simtime * self.maxField
        #const#return 0.9 * H_analytic/Constants.mu0

field = Field(maxField = 2./Constants.mu0 , simtime = 100e-9, startField = 0.9 * H_analytic/Constants.mu0)
print("Start", field.simtime, " [ns] run")
stream = open(args.outdir + "m.dat", "w")
timer = time.time()
i = 0
nvti = 0
while (state.t < field.simtime and state.mean_m(wire_dir_val) < (1. - 1e-6)):
    if i%300 == 0:
        state.write_vti(args.outdir + "m_" + str(nvti))
        nvti = nvti + 1
    if wire_dir == "z":
        fields[0].set_homogeneous_field(0.0, 0.0, field.from_time(state.t))
    elif wire_dir == "y":
        fields[0].set_homogeneous_field(0.0,  field.from_time(state.t), 0.0)
    elif wire_dir == "x":
        fields[0].set_homogeneous_field(field.from_time(state.t), 0.0, 0.0)
    if llg_before_mini == True:
        llg.step(state)
    else:
        minimizer.minimize(state)

    printzee = af.mean(af.mean(af.mean(fields[0].H_in_Apm(state), dim=0), dim=1), dim=2)

    if i%20 == 0:
        print("H[%analytic]", "{:6.4f}".format(printzee[0, 0, 0, wire_dir_val].scalar()*Constants.mu0 / H_analytic), "time[%]=", "{:6.4f}".format(state.t/field.simtime), "mean_m=", "{:5.4f}".format(state.mean_m(wire_dir_val)), "field[Oe]=", "{:6.4f}".format(field.from_time(state.t)), "field[T]=", "{:6.4f}".format(printzee[0, 0, 0, wire_dir_val].scalar()*Constants.mu0))
    stream.write("%e, %e, %e, %e, %e, %e\n" %(state.t, state.mean_m(0), state.mean_m(1), state.mean_m(2), field.from_time(state.t), printzee[0, 0, 0, wire_dir_val].scalar()))
    stream.flush()
    i = i + 1

stream.close()
state.write_vti(args.outdir + "m_switched")
print("simulated", state.t, "[s] (out of", field.simtime, ") in", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]\n")
print("H_analytic=", H_analytic, " [T]")
print("switched at Hext=", printzee[0, 0, 0, wire_dir_val].scalar()*Constants.mu0, " [T] with state.t=", state.t, " [s]")
