import sys
import arrayfire as af
from magnumaf import *
import numpy as np
import time
from math import atan, sin, cos

def get_norm_mz(a):
    x = a[0,0,0,0].scalar()
    y = a[0,0,0,1].scalar()
    z = a[0,0,0,2].scalar()
    norm = np.sqrt(x**2 + y**2 + z**2)
    return z/norm

# Setting filepath
print ("The arguments are: " , str(sys.argv))
if sys.argv[1][-1] != "/":
    sys.argv[1] = sys.argv[1] + "/"
    print ("Info: Added '/' to sys.arvg[1]: ", sys.argv[1])

filepath = sys.argv[1]

# dimensions
x, y, z = 6000e-6, 6000e-6, 2500e-6 + 1000e-6
nx, ny, nz = 100, 100, 7
dx, dy, dz = x/nx, y/ny, z/nz

# material param
Ms = 1.6/Constants.mu0

print(x, y, z, dx, dy, dz)
mesh = Mesh(nx, ny, nz, dx, dy, dz)
m = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.f64)
print(m[:, :, 0:4, :].dims(), Util.disk(nx, ny, 5, [0, 1, 0]).dims())
m[:, :, 0:5, :] = Util.disk(nx, ny, 5, [0, 1, 0])
state = State(mesh, Ms, m)
state.write_vti(filepath + "m_init")

# Initializing interaction terms
demag = DemagField(mesh, caching = True, verbose = False)
llg = LLGIntegrator(0, [demag]) # dummy object for heff readout
h = llg.h(state)
Util.write_vti(h, dx, dy, dz, filepath + "h")

# getting sensor positions on grid
ix = int((3e-3 - 1.41e-3)/(x/2) * nx/2)
iy = int((3e-3 - 1.41e-3)/(y/2) * ny/2)
iz = nz - 1 # array index start from 0
print(ix, iy, iz)

# run phi sweep
stream = open(filepath + "h.dat", "w")
ni = 90
for i in range(ni):
    phi = 2 * np.pi * (i/ni)
    state.m_partial[:, :, 0:5, :] = Util.disk(nx, ny, 5, [sin((2*np.pi) * i/ni), cos((2*np.pi) * i/ni), 0])
    state.write_vti(filepath + "m_init" + str(i))
    start = time.time()
    h = llg.h(state)
    print ("in loop", i, "with magnetization in", sin(i/ni), cos(i/ni), 0, "calc h [s]= ", time.time() - start)

    h1 = h[ ix, 0, iz, 2].scalar()
    h2 = h[ 0, iy, iz, 2].scalar()
    h3 = h[-ix, 0, iz, 2].scalar()
    h4 = h[ 0,-iy, iz, 2].scalar()

    dhx = h1 - h3
    dhy = h2 - h4

    h_phi = atan( dhy /  dhx)
    h_err = atan( dhy /  dhx) - phi
    h_err = 360/(2 * np.pi) * h_err

    m1 = get_norm_mz(h[ ix, 0, iz, 2])
    m2 = get_norm_mz(h[ 0, iy, iz, 2])
    m3 = get_norm_mz(h[-ix, 0, iz, 2])
    m4 = get_norm_mz(h[ 0,-iy, iz, 2])

    dmx = m1 - m3
    dmy = m2 - m4

    m_phi = atan( dmy /dmx )
    m_err = atan( dmy /dmx ) - phi
    m_err = 360/(2 * np.pi) * m_err


    print("h vals:", h1, h2, h3, h4, dhx, dhy, h_phi, m_phi)
    stream.write("%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n" %(phi, h_err, h_phi, m_err, m_phi, h_err, h1, h2, h3, h4, dhx, dhy) )
    stream.flush()
    Util.write_vti(h, dx, dy, dz, filepath + "h" + str(i))

stream.close()
