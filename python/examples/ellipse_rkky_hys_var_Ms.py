# Example demonstrating the MuMAG Standard Problem 4
# Run with 'magnum.af sp4.py' or 'python3 sp4.py $PWD'

import numpy as np
import arrayfire as af
from magnumaf import *
import sys
import time

af.set_device(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
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
#state = State(mesh, Ms = Ms_bottom, m = m_initi(nx, ny, nz))
state = State(mesh, Ms = Ms_initi(nx, ny, nz), m = m_initi(nx, ny, nz))
#state = State(mesh, Ms = Ms_initi(nx, ny, nz), m = m_random_sphere(nx, ny, nz))
#state = State(mesh, Ms, m = m_initi(nx, ny, nz))
state.write_vti(sys.argv[1] + "m_init")


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
stream = open(sys.argv[1] + "m.dat", "w", buffering = 1)
stream.write("# Hext [T], mx, my, mz")
for i in range(0, hys_steps + 1):
    extfield = hysteresis_factor(i, hys_steps) * ext_max_field
    ext.set_homogeneous_field(0.01 * extfield, extfield, 0) # easy-axis loop
    #hard axis loop# ext.set_homogeneous_field(extfield, 0, 0)
    if minimize:
        minimizer.minimize(state)
    else:
        llg.relax(state, precision = 1e-11, verbose = True)
    state.write_vti(sys.argv[1] + "m_step_"+ str(i))
    mx, my, mz = state.m_mean()
    print(i, 'ext[T]={:2.3f}, mx={:1.3f}, my={:1.3f}, mz={:1.3f}'.format(ext.h(state)[0, 0, 0, 0].scalar() * Constants.mu0, mx, my, mz))
    stream.write("%e, %e, %e, %e\n" %(extfield * Constants.mu0, mx, my, mz))

stream.close()

## Initial magnetization configuration
#m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
#m0[:, :, 0, 1] = -1.
#print(geom.type())
#print(RKKYarr.type())
#print(excharr.type())
#print(af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64).type())
#print(af.mean(RKKYarr))
#print(af.max(RKKYarr))
#print(af.min(RKKYarr))
#print(af.max(excharr))
#print(af.min(excharr))
#print(af.min(excharr))
# when looping over state.m use .scalar() for if
#state = State(mesh, Ms, m = Util.disk(nx, ny, nz, axis=[0, 1, 0]))
#ntrue = 0
#for i in range(0, nx):
#    for j in range(0, ny):
#        print(i, j, state.m[i, j, 1, 1].scalar())
#        if state.m[i, j, 1, 1].scalar() == 1.:
#            ntrue = ntrue + 1
#            state.m_partial[i, j, 1, 1] = af.constant(-1, 1, dtype=af.Dtype.f64)
#
#print('ntrue', ntrue)


#rkky_values = af.constant(RKKY, nx, ny, nz, 3, dtype=af.Dtype.f64)
#exch_values = af.constant(A, nx, ny, nz, 3, dtype=af.Dtype.f64)
#rkky_indices = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.u32)

## relax
#if minimize:
#    minimizer.minimize(state)
#else:
#    llg.relax(state)

#E = []
#print("Start rotating")
#stream = open(sys.argv[1]+"m.dat", "w")
#timer = time.time()
#for i in range(0, 360):
#    mix = np.cos(i * np.pi/180.);
#    miy = np.sin(i * np.pi/180.);
#    m = state.m
#    m[:, :, 1, 0] = mix
#    m[:, :, 1, 1] = miy
#    state.m = m
#    E.append( llg.E(state) )
#    print("angle=[°]", i, " E=", E[-1])
#    mean = af.mean(af.mean(af.mean(state.m, dim=0), dim=1), dim=2)
#    stream.write("%d, %e, %e, %e, %e\n" %(i, E[-1], mean[0, 0, 0, 0].scalar(), mean[0, 0, 0, 1].scalar(), mean[0, 0, 0, 2].scalar()))
#print("fullrotation in ", time.time() - timer, "[s]")
#stream.close()
#
#print("E diff =", max(E) - min(E), "[J], should be around 2.1e-19 (from mumax3 plot)")
#
## plotting data with gnuplot
#from os import system
#system('gnuplot -e "\
#    set terminal pdf;\
#    set output \'' + sys.argv[1] + 'm.pdf\';\
#    set xlabel \'angle [°]\';\
#    set ylabel \'E\';\
#    p \'' + sys.argv[1] + '/m.dat\' u 1:2 w l t \'E\';\
#    set ylabel \'<m>\';\
#    p \'' + sys.argv[1] + '/m.dat\' u 1:3 w l t \'<m_x>\',\
#    \'\' u 1:4 w l t \'<m_y>\',\
#    \'\' u 1:5 w l t \'<m_z>\';\
#"')
#
## show pdf with evince
#system('evince ' + sys.argv[1] +'m.pdf')
