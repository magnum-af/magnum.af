#!/usr/bin/python3
# Example demonstrating the MuMAG Standard Problem 4
# Run with 'magnum.af sp4.py' or 'python3 sp4.py $PWD'

import numpy as np
import arrayfire as af
from magnumaf import *
import sys
import time

args = parse()
af.info()

minimize = False # True uses minimizer, False uses llg integration
hys_steps = 120 # hysteresis steps
ext_max_field = 120e-3/Constants.mu0
# Discretization
nx, ny, nz = 100, 600, 2
x, y, z = 1e-6, 6e-6, 6e-9
dx, dy, dz = x/nx, y/ny, z/nz
print("dx={}, dy={}, dz={}".format(dx, dy, dz))

# Material parameters
Ms_bottom = 1.7/Constants.mu0
Ms_top = 1.1 * Ms_bottom
A = 13e-12
RKKY = -0.1e-3 * dz

def m_random_sphere(nx, ny, nz):
    m = np.zeros((nx, ny, nz, 3))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(r<=1):
                for iz in range(0, nz):
                    randx = np.random.normal()
                    randy = np.random.normal()
                    randz = np.random.normal()
                    #print(randx, randy, randz)
                    norm = sqrt(randx**2+randy**2+randz**2)
                    randx = randx/norm
                    randy = randy/norm
                    randz = randz/norm
                    m[ix, iy, iz, 0] = randx
                    m[ix, iy, iz, 1] = randy
                    m[ix, iy, iz, 2] = randz
    #print(np.mean(np.mean(np.mean( m, axis=0) , axis=1) , axis=2))
    #print(np.mean(np.mean(np.mean(np.mean( m, axis=0) , axis=1) , axis=2) , axis=3))
    return af.from_ndarray(m)

def m_initi(nx, ny, nz):
    m = np.zeros((nx, ny, nz, 3))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(r<=1):
                m[ix, iy, 0, 1] =  1.
                m[ix, iy, 1, 1] = -1.
                # tiling m[ix, iy, 0, 0] = 0.001
                # tiling m[ix, iy, 1, 0] = 0.001
    return af.from_ndarray(m)

def Ms_initi(nx, ny, nz):
    Ms = np.zeros((nx, ny, nz))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(r<=1):
                Ms[ix, iy, 0] = Ms_bottom
                Ms[ix, iy, 1] = Ms_top
    return af.from_ndarray(Ms)

def disk(nx, ny, nz):
    m = np.zeros((nx, ny, nz))
    for ix in range (0, nx):
        for iy in range(0, ny):
            a= nx/2
            b= ny/2
            rx=ix-nx/2.
            ry=iy-ny/2.
            r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
            if(r<=1):
                m[ix, iy, :] = 1.
    return af.from_ndarray(m)

geom = disk(nx, ny, nz)
RKKYarr = (geom == 1) * RKKY
excharr = (geom == 1) * A

# Creating objects
mesh = Mesh(nx, ny, nz, dx, dy, dz)
state = State(mesh, Ms = Ms_initi(nx, ny, nz), m = m_initi(nx, ny, nz))
state.write_vti(args.dir + "m_init")


demag = DemagField(mesh, verbose = True, caching = True, nthreads = 6)
rkkyexch = RKKYExchangeField(RKKYarr, excharr, mesh, rkky_indices = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.u32))


ext_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
ext = ExternalField(ext_field)
terms = [demag, rkkyexch, ext]

if minimize:
    minimizer = LBFGS_Minimizer(terms, tol=1e-15, maxiter=1000)
else:
    llg = LLGIntegrator(alpha = 1, terms = terms)

def hysteresis_factor(i, steps):
    if i < steps/4.:
        return 4 * i/steps
    elif i < 3 * steps/4.:
        return (4 * -i/steps + 2. )
    else:
        return (4 * i/steps - 4. )


# running hysteresis loop
stream = open(args.dir + "m.dat", "w", buffering = 1)
stream.write("# Hext [T], mx, my, mz")
for i in range(0, hys_steps + 1):
    extfield = hysteresis_factor(i, hys_steps) * ext_max_field
    ext.set_homogeneous_field(0.01 * extfield, extfield, 0) # easy-axis loop
    #hard axis loop# ext.set_homogeneous_field(extfield, 0, 0)
    if minimize:
        minimizer.minimize(state)
    else:
        llg.relax(state, precision = 1e-11, verbose = True)
    state.write_vti(args.dir + "m_step_"+ str(i))
    mx, my, mz = state.mean_m()
    print(i, 'ext[T]={:2.3f}, mx={:1.3f}, my={:1.3f}, mz={:1.3f}'.format(ext.H_in_Apm(state)[0, 0, 0, 0].scalar() * Constants.mu0, mx, my, mz))
    stream.write("%e, %e, %e, %e\n" %(extfield * Constants.mu0, mx, my, mz))

stream.close()
