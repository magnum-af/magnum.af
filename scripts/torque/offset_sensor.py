# argv[1] filepath
# argv[2] GPU number [0-3] (0)
# argv[3] simtime [ns] (1)

import arrayfire as af
from magnumaf import *
import sys
import time

if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))
af.info()

# Physical dimensions in [m]
x = 1e-9
y = 1e-9
z = 0.6e-9
# Discretization
nx = 2
ny = 2
nz = 2

mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
material = Material(ms = 1.4e6, A = 30e-12, alpha = 0.02, Ku1=1e4, Ku1_axis=[1,0,0])
m = Util.normed_homogeneous_field(nx, ny, nz, [1,0,0])
state = State(mesh, material, m)

polarization = Util.normed_homogeneous_field(nx, ny, nz, [0, 1, 0])

fields = [
    DemagField(mesh, material),
    ExchangeField(mesh, material),
    SpinTransferTorqueField(polarization, nu_damp=.3, nu_field=.4, j_e=2e10),
    UniaxialAnisotropyField(mesh, material),
    ExternalField(Util.normed_homogeneous_field(nx, ny, nz, [1,1,0], 10e-3/Constants.mu0)),
]
print (fields)
Llg = LLGIntegrator(terms=fields, mode="RKF45", hmin = 1e-15, hmax = 3.5e-10, atol = 1e-6, rtol = 1e-6)

#print(fields[2].polarization_field)
#print(state.m)

simtime = float(sys.argv[3])*1e-9 if len(sys.argv) > 3 else 1e-9
nstep = int(sys.argv[4]) if len(sys.argv) > 4 else 100

stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
itcount=0
print ("Starting integration for", simtime, " [ns]. Calculating mean every ", nstep, "step.")
while state.t < simtime:
    Llg.step(state)
    if itcount % nstep == 0:
        mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2) # calculates spacially averaged <mx>, <my>, <mz>
        print(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar())
        stream.write("%e, %e, %e, %e\n" %(state.t, mean[0,0,0,0].scalar(), mean[0,0,0,1].scalar(), mean[0,0,0,2].scalar()))
        stream.flush()
    itcount=itcount+1
print("Simulated ", simtime, " [ns] in ", time.time() - timer, "[s]")
stream.close()


#heff = llg.h(state)
#print(heff[0,0,0,0].scalar())
#print(heff[0,0,0,1].scalar())
#print(heff[0,0,0,2].scalar())
