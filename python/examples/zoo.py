import arrayfire as af
import numpy as np
from magnumaf import *
import sys
import time


af.set_device(Util.gto_gpu_renumeration(int(sys.argv[2])) if len(sys.argv) > 2 else 0)
af.info()

# Physical dimensions in [m]
x, y, z =  3e-6, 3e-6, 62e-9
# Discretization
nx, ny, nz = 100, 100, 4

A = 1e-11
K = 1e4
Ms = 320e3
D_max = 2.5e-3 # DMI max
H_max = 0.2 # H max in [T]

print('sqrt(A/K)=', np.sqrt(A/K))

dn = 10
start = time.time()
for iD in range(0, dn + 1):
    for iH in range(0, dn + 1):
        timer = time.time()
        D = (iD/dn) * D_max
        H = (iH/dn) * H_max
        print('starting step: D=', D, ', H=', H)

        # Creating objects
        mesh = Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)
        
        # Initial magnetization configuration
        m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
        m0[:, :, :, 2] = 1.
        
        state = State(mesh, Ms = Ms, m = m0)
        demag = DemagField(mesh, verbose = False, caching = True, nthreads = 6)
        exch = ExchangeField(A)
        dmi = DmiField(D, [0, 0, 1])
        
        zeeswitch = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
        zeeswitch[:, :, :, 2] = H/Constants.mu0
        zee = ExternalField(zeeswitch)
        
        llg = LLGIntegrator(alpha = 1, terms = [demag, exch, dmi, zee])
        
        # Relaxing
        llg.relax(state, precision = 5e-8, ncalcE = 100, nprint = 1000, verbose = True)
        state.write_vti(sys.argv[1] + 'm_relaxed_iD' + str(iD) + '_iH' + str(iH))
        print('step', (iD * (dn + 1) + iH), '/', (dn + 1)**2 - 1 , " took ", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]")
