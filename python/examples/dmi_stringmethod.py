# Example demonstrating the MuMAG Standard Problem 4
# Run with 'magnum.af sp4.py' or 'python3 sp4.py $PWD'

import numpy as np
import arrayfire as af
from magnumaf import *
import sys
import time

af.info()
start = time.time()

# micromagnetic parameters
Ms = 1.1e6;
A = 1.6e-11;
D = 2 * 5.76e-3;
Ku1 = 6.4e6;
extfield = 10e-3 / Constants.mu0;

# Discretization
nx, ny, nz = 30, 30, 1
dx = 1e-9

# atomistic parameters
J_atom = 2. * A * dx;
D_atom = D * pow(dx, 2);
K_atom = Ku1 * pow(dx, 3);
p = Ms * pow(dx, 3);


# Creating objects
mesh = Mesh(nx, ny, nz, dx, dx, dx)

def disk(nx, ny, nz, r = 0.5):
    m = np.zeros((nx, ny, nz, 3))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            rr = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(rr <= r):
                m[ix, iy, :, :] =-1.
            else:
                m[ix, iy, :, :] = 1.
    return af.from_ndarray(m)

state1 = State(mesh, Ms = p, m = disk(nx, ny, nz))
state1.write_vti(sys.argv[1] + "m_init")

zee = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
zee[:, :, :, 2] = extfield

exch = AtomisticExchangeField(J_atom)
dmi  = AtomisticDmiField(D_atom, [0, 0, -1])
ani  = AtomisticUniaxialAnisotropyField(K_atom, [0, 0, 1])
dip  = AtomisticDipoleDipoleField(mesh)
ext  = AtomisticExternalField(zee)
terms = [exch, dmi, ani, dip, ext]

llg = LLGIntegrator(alpha = 1, terms = terms)

print("relaxing 1ns")
stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
while state1.t < 1e-09:
    llg.step(state1)
    mean = af.mean(af.mean(af.mean(state1.m, dim=0), dim=1), dim=2)
    stream.write("%e, %e, %e, %e\n" %(state1.t, mean[0, 0, 0, 0].scalar(), mean[0, 0, 0, 1].scalar(), mean[0, 0, 0, 2].scalar()))
stream.close()
print("relaxed in", time.time() - timer, "[s]")
state1.write_vti(sys.argv[1] + "m_relaxed")

# Initial magnetization configuration
m2 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m2[:, :, :, 2] = 1
state2 = State(mesh, Ms = p, m = m2)
state2.write_vti(sys.argv[1] + "m_init")
print("relaxing 1ns")
stream = open(sys.argv[1]+"m2.dat", "w")
timer = time.time()
while state2.t < 1e-09:
    llg.step(state2)
    mean = af.mean(af.mean(af.mean(state2.m, dim=0), dim=1), dim=2)
    stream.write("%e, %e, %e, %e\n" %(state2.t, mean[0, 0, 0, 0].scalar(), mean[0, 0, 0, 1].scalar(), mean[0, 0, 0, 2].scalar()))
stream.close()
print("relaxed in", time.time() - timer, "[s]")
state2.write_vti(sys.argv[1] + "m_relaxed")

n_interp = 60
string_dt = 5e-14

string = String(state1, [state1.m, state2.m], n_interp, string_dt, llg)
# also works #string = String(state1, [state1, state2], n_interp, string_dt, llg)
print("initialized string")
dE = string.run(sys.argv[1], verbose = False)
expected_result = 1.045368e-19
print("dE=", dE, "expected", expected_result)

print("total time =", time.time() - start, "[s]")


# plotting data with gnuplot
from os import system
system('gnuplot -e "\
    set terminal pdf;\
    set output \'' + sys.argv[1] + 'm.pdf\';\
    set xlabel \'t [ns]\';\
    set ylabel \'<m>\';\
    p \'' + sys.argv[1] + '/m.dat\' u 1:2 w l t \'<m_x>\',\
    \'\' u 1:3 w l t \'<m_y>\',\
    \'\' u 1:4 w l t \'<m_z>\';\
"')

# show pdf with evince
system('evince ' + sys.argv[1] +'m.pdf')
