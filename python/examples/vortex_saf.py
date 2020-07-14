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
x, y, z =  1000e-9, 1000e-9, 84e-9
# Discretization
nx, ny, nz = 250, 250, 42
dx, dy, dz =x/nx, y/ny, z/nz

Ms = 1.393e6 # //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
A = 1.5e-11 #  //[J/m]
ext_in_bottom_ref = 275e-3/Constants.mu0 #

# when defining E = - 2 RKKY \int_{A} m_a m_b dA
# other definitions skip the 2 s.t. RKKY_alternative = 2 * RKKY
RKKY_val = -0.4e-3 * dz # [J/m^2] * [m]

# Creating objects
mesh = Mesh(nx, ny, nz, dx=x/nx, dy=y/ny, dz=z/nz)
print(mesh.nx, mesh.ny, mesh.nz, mesh.dx, mesh.dy, mesh.dz, dx, dy, dz)

# Initial magnetization configuration
nvort = 39 # number of layer in vortex
diffnz = nz - nvort

timer = time.time()
vortex = Util.vortex(Mesh(nx, ny, nz-diffnz, dx=x/nx, dy=y/ny, dz=z/nz))
print("Init vortex in", time.time() - timer, "[s]")# TODO use cython to speedup
#print(vortex.dims())

Util.write_vti(vortex, dx, dy, dz, filepath + "vortex_init")
m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, 0, 0] = 1.
m0[:, :, 1, 0] = -1.
m0[:, :, diffnz:, :] = vortex

#Util.write_vti(m0, dx, dy, dz, filepath + "m0")

m0 = Util.normalize(m0)

demag = DemagField(mesh, verbose = True, caching = True, nthreads = 16)

def disk(nx, ny, nz):
    m = np.zeros((nx, ny, nz, 1))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(r<=1):
                m[ix, iy, :, :] = 1.
    return af.from_ndarray(m)

geom = af.tile(disk(nx, ny, nz), 1, 1, 1, 3)
geom[:, :, 0, :] = 1.
geom[:, :, 1, :] = 1.
geom[:, :, 2, :] = 0.
Util.write_vti(geom, dx, dy, dz, filepath + "geom")
excharr = (geom == 1) * A

RKKYarr = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
RKKYarr[:, :, 0, :] = RKKY_val
RKKYarr[:, :, 1, :] = RKKY_val
Util.write_vti(excharr, dx, dy, dz, filepath + "excharr")
Util.write_vti(RKKYarr, dx, dy, dz, filepath + "RKKYarr")
rkkyexch = RKKYExchangeField(RKKYarr, excharr, mesh, rkky_indices = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.u32))

region_ref = disk(nx, ny, 1) # region to evaluate mean values in reference layer
Util.write_vti(region_ref, dx, dy, dz, filepath + "region_ref")

stream_angle = open(filepath + "angle_over_h.dat", "w")
stream_angle.write("#ext_y[mT]\ttile along x-axis[°]\n")

for iext_y in range(0,6):
    ext_y= 20 * iext_y # field in mT
    filepath=filepath_const + '/' + str(ext_y) + 'mT/'
    try:
        os.mkdir(filepath)
    except:
        print("mkdir failed")
    print(filepath)
    state = State(mesh, Ms, m = m0)
    print('state.mean_m()', state.mean_m())
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
    #while state.t < 5e-12:
    while state.t < 5e-9:
        llg.step(state)
        mx, my, mz = state.mean_m()
        print('step={:d}, t[ns]={:1.6e}, mx={:1.6f}, my={:1.6f}, mz={:1.6f}'.format(llg.accumulated_steps, state.t * 1e9, mx, my, mz))
        stream.write('{:d}\t{:1.6e}\t{:1.6f}\t{:1.6f}\t{:1.6f}\n'.format(llg.accumulated_steps, state.t * 1e9, mx, my, mz))
        if llg.accumulated_steps % write_vti_every == 0:
            state.write_vti(filepath + "m_relaxing_step_" + str(llg.accumulated_steps))
            Util.write_vti(demag.h(state), dx, dy, dz, filepath + "demag_step_" + str(llg.accumulated_steps))
    print("relaxed in", time.time() - timer, "[s]")
    state.write_vti(filepath + "m_relaxed")

    af_m_layer0 = state.m[:, :, 0, :]
    af_m_layer1 = state.m[:, :, 1, :]
    af_m_layer3 = state.m[:, :, 3, :]
    np_m_layer0 = af_m_layer0.to_ndarray()
    np_m_layer1 = af_m_layer1.to_ndarray()
    np_m_layer3 = af_m_layer3.to_ndarray()
    #print(np_m_layer0.shape)
    #print(np_m_layer1.shape)
    #print(np_m_layer3.shape)
    np.savetxt(filepath + 'layer0_0.txt', np_m_layer0[:, :, 0, 0])
    np.savetxt(filepath + 'layer0_1.txt', np_m_layer0[:, :, 0, 1])
    np.savetxt(filepath + 'layer0_2.txt', np_m_layer0[:, :, 0, 2])

    np.savetxt(filepath + 'layer1_0.txt', np_m_layer1[:, :, 0, 0])
    np.savetxt(filepath + 'layer1_1.txt', np_m_layer1[:, :, 0, 1])
    np.savetxt(filepath + 'layer1_2.txt', np_m_layer1[:, :, 0, 2])

    np.savetxt(filepath + 'layer3_0.txt', np_m_layer3[:, :, 0, 0])
    np.savetxt(filepath + 'layer3_1.txt', np_m_layer3[:, :, 0, 1])
    np.savetxt(filepath + 'layer3_2.txt', np_m_layer3[:, :, 0, 2])

    pinned_layer = state.m[:, :, 0, :]
    Util.write_vti(pinned_layer, dx, dy, dz, filepath + "pinned_layer")
    reference_layer = state.m[:, :, 1, :]
    Util.write_vti(reference_layer, dx, dy, dz, filepath + "reference_layer")
    mean_ref = Util.spacial_mean_in_region(reference_layer, region_ref) # evaluating mean_m only in circlic region under vortex
    print('mean pinned       =', af.mean(af.mean(pinned_layer, dim=0), dim=1).to_ndarray()[0, 0, 0, :])
    print('mean ref          =', af.mean(af.mean(reference_layer, dim=0), dim=1).to_ndarray()[0, 0, 0, :])
    print('mean ref in region=', mean_ref)
    tile_in_degree = 360/(2 * np.pi) * np.arctan(mean_ref[1]/mean_ref[0])
    print('arctan[°]=', tile_in_degree)
    stream_angle.write('{:1.6f}\t{:1.6f}\n'.format(ext_y, tile_in_degree))

stream.close()
stream_angle.close()
