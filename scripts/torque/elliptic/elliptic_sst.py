# argv[1] filepath
# argv[2] GPU number [0-3] (0)

import arrayfire as af
from magnum_af import *
import sys
import time

if len(sys.argv) > 2:
    af.set_device(int(sys.argv[2]))
af.info()

# Physical dimensions in [m]
x = 700e-9
y = 4000e-9
z = 1.6e-9
# Discretization 
nx = int(x/5e-9)
ny = int(y/5e-9)
nz = 1

print(nx, ny, nz, x, y, z)

m, n_cells = Util.disk(nx, ny, nz, [1,1,0])
material = Material(ms = 8.6e5, A = 30e-12, alpha = 0.1)
mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz)
state = State(mesh, material, m)

polarization = Util.normed_homogeneous_field(nx, ny, nz, [1, 0, 0]) # Current in pinned layer along y-axis creates polarization in ellipse (which is in positive z dir) in (+/-)? x-dir

fields = [
    DemagField(mesh, material, verbose = True, caching = True),
    ExchangeField(mesh, material),
    SpinTransferTorqueField(polarization, nu_damp=.1, nu_field=.7, j_e=1.6e11),
    #UniaxialAnisotropyField(mesh, material),
    #Zee(Util.normed_homogeneous_field(nx, ny, nz, [1,1,0], 10e-3/Constants.mu0)),
]
print (fields)

# LLG version
llg = LLGIntegrator(terms=fields, mode="RKF45", hmin = 1e-15, hmax = 3.5e-10, atol = 1e-6, rtol = 1e-6)
#llg.relax(state)

stream = open(sys.argv[1]+"m.dat", "w")
timer = time.time()
itcount=0
simtime = float(sys.argv[3])*1e-9 if len(sys.argv) > 3 else 1e-9
nstep = int(sys.argv[4]) if len(sys.argv) > 4 else 100
print ("Starting integration for", simtime, " [ns]. Calculating mean every ", nstep, "step.")
while state.t < simtime:
    llg.llgstep(state)
    print(state.t, state.meanxyz(0), state.meanxyz(1), state.meanxyz(2), state.steps)
    stream.write("%e, %e, %e, %e, %d\n" %(state.t, state.meanxyz(0), state.meanxyz(1), state.meanxyz(2), state.steps))
    stream.flush()
    if itcount % nstep == 0:
        state.py_vti_writer_micro(sys.argv[1] + "step" + str(itcount))
    itcount=itcount+1
print("Simulated ", simtime, " [ns] in ", time.time() - timer, "[s]")
stream.close()

# Minimizer version
#timer = time.time()
#minimizer = LBFGS_Minimizer(terms=fields, tol=1e-15, maxiter=1000)
#minimizer.pyMinimize(state)
#print("Minimized in ", time.time() - timer, "[s]")
#state.py_vti_writer_micro(sys.argv[1] + "minimized")

mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
print("Mean magnetization: ",  state.meanxyz(0), state.meanxyz(1), state.meanxyz(2))
