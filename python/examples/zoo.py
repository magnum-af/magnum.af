#!/usr/bin/python3
import arrayfire as af
import numpy as np
from magnumaf import *
import sys
import time

args = parse()
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

        #llg.relax(state, precision = 5e-8, ncalcE = 100, nprint = 1000, verbose = True)

        E_prev = 1e20;
        precision = 1e-8
        integration_time = 5e-9
        steps_per_iteration = 100
        print_steps_every = 10
        write_vti_every = 100
        while state.t < integration_time:
            E_curr = llg.Eeff_in_J(state)
            Ediff = abs((E_prev - E_curr) / E_prev)
            print('Ediff', Ediff)
            if Ediff < precision:
                break
            E_prev = E_curr
        
            for i in range(steps_per_iteration):
                llg.step(state)
                mx, my, mz = state.mean_m()
                if llg.accumulated_steps % print_steps_every == 0:
                    print('step={:d}, t[ns]={:1.6e}, mx={:1.6f}, my={:1.6f}, mz={:1.6f}'.format(llg.accumulated_steps, state.t * 1e9, mx, my, mz))
                if llg.accumulated_steps % write_vti_every == 0:
                    state.write_vti(args.outdir + "m_iD" + str(iD) + '_iH' + str(iH) + "step_" + str(llg.accumulated_steps) )

        # Relaxing
        state.write_vti(args.outdir + 'm_relaxed_iD' + str(iD) + '_iH' + str(iH))
        print('step', (iD * (dn + 1) + iH), '/', (dn + 1)**2 - 1 , " took ", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]")
