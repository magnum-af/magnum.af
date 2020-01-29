# Example demonstrating the MuMAG Standard Problem 4
# Run with 'magnum.af sp4.py' or 'python3 sp4.py $PWD'

import os
import numpy as np
import arrayfire as af
from magnumaf import *
import sys
import time

af.set_device(Util.gto_gpu_renumeration(int(sys.argv[2])) if len(sys.argv) > 2 else 0)
af.info()
start = time.time()
filepath_const=sys.argv[1]
filepath=filepath_const

# Physical dimensions in [m]
x, y, z =  1000e-9, 1000e-9, 4e-9
# Discretization
nx, ny, nz = 250, 250, 2
dx, dy, dz =x/nx, y/ny, z/nz

Ms = 1.393e6 # //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
A = 1.5e-11 #  //[J/m]
ext_in_bottom_ref = 300e-3/Constants.mu0 #
RKKY_val = -0.4e-3 * dz

# Creating objects
mesh = Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)
print(mesh.nx, mesh.ny, mesh.nz, mesh.dx, mesh.dy, mesh.dz, dx, dy, dz)

m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, 0, 0] = 1.
m0[:, :, 1, 0] = -1.

demag = DemagField(mesh, verbose = True, caching = True, nthreads = 16)

excharr = af.constant(A, nx, ny, nz, 3, dtype=af.Dtype.f64)
RKKYarr = af.constant(RKKY_val, nx, ny, nz, 3, dtype=af.Dtype.f64)
Util.write_vti(excharr, dx, dy, dz, filepath + "excharr")
Util.write_vti(RKKYarr, dx, dy, dz, filepath + "RKKYarr")
rkkyexch = RKKYExchangeField(RKKYarr, excharr, mesh, rkky_indices = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.u32))

for iext_y in range(0,6):
    ext_y= 20 * iext_y # field in mT
    filepath=filepath_const + '/' + str(ext_y) + 'mT/'
    try:
        os.mkdir(filepath)
    except:
        print("mkdir failed")
    print(filepath)
    state = State(mesh, Ms, m = m0)
    print('state.m_mean()', state.m_mean())
    state.write_vti(filepath + "state_m0")
    #Util.write_vti(state.Ms, dx, dy, dz, filepath + "Ms")

    extf = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
    extf[:, :, 0, 0] = ext_in_bottom_ref

    ext_global = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
    print('(ext_y * 1e-3)/Constants.mu0=', (ext_y * 1e-3)/Constants.mu0)
    ext_global[:, :, :, 1] = (ext_y * 1e-3)/Constants.mu0
    ext = ExternalField(extf + ext_global)
    Util.write_vti(extf + ext_global, dx, dy, dz, filepath + "ext_combined")

    terms = [demag, rkkyexch, ext]
    llg = LLGIntegrator(alpha = 1, terms = terms)

    # Relaxing
    print("relaxing 1ns")
    stream = open(filepath + "m.dat", "w")
    stream.write('#step={:d}, t[ns]={:1.6e}, mx={:1.6f}, my={:1.6f}, mz={:1.6f}\n')
    timer = time.time()
    write_vti_every = 100
    while state.t < 1e-9:
        llg.step(state)
        mx, my, mz = state.m_mean()
        #print('step={:d}, t[ns]={:1.6e}, mx={:1.6f}, my={:1.6f}, mz={:1.6f}'.format(llg.accumulated_steps, state.t * 1e9, mx, my, mz))
        stream.write('{:d}\t{:1.6e}\t{:1.6f}\t{:1.6f}\t{:1.6f}\n'.format(llg.accumulated_steps, state.t * 1e9, mx, my, mz))
        if llg.accumulated_steps % write_vti_every == 0:
            state.write_vti(filepath + "m_relaxing_step_" + str(llg.accumulated_steps))
            Util.write_vti(demag.h(state), dx, dy, dz, filepath + "demag_step_" + str(llg.accumulated_steps))
    print("relaxed in", time.time() - timer, "[s]")
    state.write_vti(filepath + "m_relaxed")

    af_m_layer0 = state.m[:, :, 0, :]
    af_m_layer1 = state.m[:, :, 1, :]
    np_m_layer0 = af_m_layer0.to_ndarray()
    np_m_layer1 = af_m_layer1.to_ndarray()
    np.savetxt(filepath + 'layer0_0.txt', np_m_layer0[:, :, 0, 0])
    np.savetxt(filepath + 'layer0_1.txt', np_m_layer0[:, :, 0, 1])
    np.savetxt(filepath + 'layer0_2.txt', np_m_layer0[:, :, 0, 2])

    np.savetxt(filepath + 'layer1_0.txt', np_m_layer1[:, :, 0, 0])
    np.savetxt(filepath + 'layer1_1.txt', np_m_layer1[:, :, 0, 1])
    np.savetxt(filepath + 'layer1_2.txt', np_m_layer1[:, :, 0, 2])

    pinned_layer = state.m[:, :, 0:2, :]
    Util.write_vti(pinned_layer, dx, dy, dz, filepath + "pinned_layer")
    meanlist=af.mean(af.mean(pinned_layer, dim=0), dim=1).to_ndarray()
    print(meanlist)
    saf_l1_x = meanlist[0, 0, 1, 0]
    saf_l1_y = meanlist[0, 0, 1, 1]
    print("H [mT]=", ext_y, ", Angle saf layer 1 [rad]=", np.arctan(saf_l1_y/saf_l1_x))
    print("H [mT]=", ext_y, ", Angle saf layer 1 [deg]=", - (360/(2 *np.pi)) * np.arctan(saf_l1_y/saf_l1_x))
    print("")
