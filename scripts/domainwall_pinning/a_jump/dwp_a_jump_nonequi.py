import arrayfire as af
from magnumaf import *
from math import sqrt
import sys
import time

if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))
af.info()
start = time.time()

# setting wire direction (x, y or z)
wire_dir = (sys.argv[3]) if len(sys.argv) > 3 else "z"
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

# Discretization

# Material
soft_Aex        = 1.0e-11/3.0
soft_ms         = 1. / Constants.mu0
soft_K_uni      = 1e5

hard_Aex        = 1.0e-11
hard_ms         = soft_ms
hard_K_uni      = soft_K_uni

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
m = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
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
z_spacing = []
length=0
factor_cell2_vs_cell1 = 1 # range{0.5, inf}
for i in range(nz):
    print("i=", i)
    if nz == 1:
        z_spacing.append(z/nz)
        length = length + z/nz
    elif i % 2 == 0:
        zval = (1/factor_cell2_vs_cell1) * z/nz
        z_spacing.append(zval)
        length = length + zval
    else:
        zval = (2 - 1/factor_cell2_vs_cell1) * z/nz
        z_spacing.append(zval)
        length = length + zval
ne_mesh = NonequispacedMesh(nx, ny, x/nx, y/ny, z_spacing)
print("length=", length)
print(ne_mesh.nz)
print(ne_mesh.z_spacing)

# Setting A, Ms and Ku1 values as fields
A_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
Ms_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)
Ku1_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)

if wire_dir == "z":
    A_field  [:, :, :nz/2, :] = soft_Aex # soft
    A_field  [:, :, nz/2:, :] = hard_Aex # hard
    Ms_field [:, :, :nz/2, :] = soft_ms # soft
    Ms_field [:, :, nz/2:, :] = hard_ms # hard
    Ku1_field[:, :, :nz/2, :] = soft_K_uni # soft
    Ku1_field[:, :, nz/2:, :] = hard_K_uni # hard
    Ku1_axis=[0, 0, 1]
elif wire_dir == "y":
    A_field  [:, :ny/2, :, :] = soft_Aex # soft
    A_field  [:, ny/2:, :, :] = hard_Aex # hard
    Ms_field [:, :ny/2, :, :] = soft_ms # soft
    Ms_field [:, ny/2:, :, :] = hard_ms # hard
    Ku1_field[:, :ny/2, :, :] = soft_K_uni # soft
    Ku1_field[:, ny/2:, :, :] = hard_K_uni # hard
    Ku1_axis=[0, 1, 0]
elif wire_dir == "x":
    A_field  [:nx/2, :, :, :] = soft_Aex # soft
    A_field  [nx/2:, :, :, :] = hard_Aex # hard
    Ms_field [:nx/2, :, :, :] = soft_ms # soft
    Ms_field [nx/2:, :, :, :] = hard_ms # hard
    Ku1_field[:nx/2, :, :, :] = soft_K_uni # soft
    Ku1_field[nx/2:, :, :, :] = hard_K_uni # hard
    Ku1_axis=[1, 0, 0]

state = State(mesh, Ms_field, m)
state.nonequimesh = ne_mesh #TODO should be handley in more object oriented way
state.write_vti(sys.argv[1] + "minit")

fields = [
    ExternalField(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f32)),
    NonequiExchangeField(A_field, ne_mesh, verbose = True),
    #SparseExchangeField(A_field, mesh),
    UniaxialAnisotropyField(Ku1_field, Ku1_axis),
]
Llg = LLGIntegrator(alpha=1.0, terms=fields)


##TODO# maxField = 0. # Note: for zero field, the domainwall does not pin to the interface for current implementation i.e. sparse yielding correct Hanalytic

class Field:
    def __init__(self, maxField = 2./Constants.mu0 , simtime = 100e-9, startField = 0/Constants.mu0):
        self.maxField = maxField # [Oe]
        self.simtime = simtime # [s]
        self.startField = startField # [Oe]
    def from_time(self, time):
        return self.startField + time / self.simtime * self.maxField

field = Field(maxField = 2./Constants.mu0 , simtime = 100e-9, startField = 0.9 * H_analytic/Constants.mu0)
print("Start", field.simtime, " [ns] run")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
i = 0
nvti = 0
while (state.t < field.simtime and state.m_mean(wire_dir_val) < (1. - 1e-6)):
    if i%300 == 0:
        state.write_vti(sys.argv[1] + "m_" + str(nvti))
        nvti = nvti + 1
    if wire_dir == "z":
        fields[0].set_homogeneous_field(0.0, 0.0, field.from_time(state.t))
    elif wire_dir == "y":
        fields[0].set_homogeneous_field(0.0,  field.from_time(state.t), 0.0)
    elif wire_dir == "x":
        fields[0].set_homogeneous_field(field.from_time(state.t), 0.0, 0.0)
    Llg.step(state)
    printzee = af.mean(af.mean(af.mean(fields[0].h(state), dim=0), dim=1), dim=2)

    if i%20 == 0:
        print("H[%analytic]", "{:6.4f}".format(printzee[0, 0, 0, wire_dir_val].scalar()*Constants.mu0 / H_analytic), "time[%]=", "{:6.4f}".format(state.t/field.simtime), "mean_m=", "{:5.4f}".format(state.m_mean(wire_dir_val)), "field[Oe]=", "{:6.4f}".format(field.from_time(state.t)), "field[T]=", "{:6.4f}".format(printzee[0, 0, 0, wire_dir_val].scalar()*Constants.mu0))
    stream.write("%e, %e, %e, %e, %e, %e\n" %(state.t, state.m_mean(0), state.m_mean(1), state.m_mean(2), field.from_time(state.t), printzee[0, 0, 0, wire_dir_val].scalar()))
    stream.flush()
    i = i + 1

stream.close()
state.write_vti(sys.argv[1] + "m_switched")
print("simulated", state.t, "[s] (out of", field.simtime, ") in", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]\n")
print("H_analytic=", H_analytic, " [T]")
print("switched at Hext=", printzee[0, 0, 0, wire_dir_val].scalar()*Constants.mu0, " [T] with state.t=", state.t, " [s]")
