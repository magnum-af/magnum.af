#  RKKY-coupled layers with comparison to analytical solution

import numpy as np
import arrayfire as af
from magnumaf import *
import sys
import time

af.set_device(Util.gto_gpu_renumeration(int(sys.argv[2])) if len(sys.argv) > 2 else 0)
af.info()

minimize = False # True uses minimizer, False uses llg integration
hys_steps = 1000 # hysteresis steps
precision = 1e-10
#hys_steps = 2000 # hysteresis steps

# Discretization
nx, ny, nz = 10, 10, 2
x, y, z = 1e-8, 1e-8, 2e-9
dx, dy, dz = x/nx, y/ny, z/nz
print("dx={}, dy={}, dz={}".format(dx, dy, dz))

#def m_analytical(dz, Ms, H, RKKY):
#    return np.sin(np.arctan( (dz * Ms * H)/(2 * abs(RKKY))))

def phi_analytical(dz, Ms, H, RKKY):
    #return np.arctan( (dz * Ms * H)/(4 * abs(RKKY)))
    return np.arctan( (dz * Ms * H)/(2 * abs(RKKY)))
    #return np.arctan( (dz * Ms * H)/(abs(RKKY)))

def mx_analytical(dz, Ms, H, RKKY):
    return np.cos(phi_analytical(dz, Ms, H, RKKY))

def my_analytical(dz, Ms, H, RKKY):
    return np.sin(phi_analytical(dz, Ms, H, RKKY))

#def H_analytical(RKKY_surface, dz, Ms):
#    return (4. * abs(RKKY_surface))/(dz * Ms)

# Material parameters
Js = 1 # [T]
Ms = Js/Constants.mu0 # "[T/mu0]"
K_pinned = 1e10
RKKY_surface = -0.5e-3
RKKY = RKKY_surface * dz
#analytical_max_field = H_analytical(RKKY_surface, dz, Ms)
#print('Analytical max field [T]  =', analytical_max_field)
#print('Analytical max field [A/m]=', analytical_max_field/Constants.mu0)
ext_max_field = 4.100/Constants.mu0
#ext_max_field = 1e0 * analytical_max_field/Constants.mu0
#ext_max_field = 1.1 * analytical_max_field/Constants.mu0

# Creating objects
mesh = Mesh(nx, ny, nz, dx, dy, dz)

m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
m0[:, :, 0, 0] =-1.
m0[:, :, 1, 0] = 1.
state = State(mesh, Ms, m = m0)
state.write_vti(sys.argv[1] + "m_init")


RKKYarr = af.constant(RKKY, nx, ny, nz, 3, dtype=af.Dtype.f64)
excharr = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.f64)
rkkyexch = RKKYExchangeField(RKKYarr, excharr, mesh, rkky_indices = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.u32))

Karr = af.constant(0, nx, ny, nz, 3, dtype=af.Dtype.f64)
Karr[:, :, 1, :] = K_pinned
aniso = UniaxialAnisotropyField(Karr, [1, 0, 0])

ext_field = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
ext = ExternalField(ext_field)
terms = [rkkyexch, aniso, ext]

if minimize:
    minimizer = LBFGS_Minimizer(terms, tol=1e-15, maxiter=1000)
else:
    llg = LLGIntegrator(alpha = 1, terms = terms, hmax = 3.5e-10, dissipation_term_only = True)

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
for i in range(1, hys_steps + 1):
    if i > hys_steps/4:
        break
    extfield = hysteresis_factor(i, hys_steps) * ext_max_field
    ext.set_homogeneous_field(0, extfield, 0)
    if minimize:
        minimizer.minimize(state)
    else:
        llg.relax(state, precision, verbose = False)
    state.write_vti(sys.argv[1] + "m_step_"+ str(i))
    mx, my, mz = state.m_mean()
    mx_top = af.mean(af.mean(state.m[:, :, 1, 0], dim=0), dim=1).scalar()
    mx_bottom = af.mean(af.mean(state.m[:, :, 0, 0], dim=0), dim=1).scalar()
    my_top = af.mean(af.mean(state.m[:, :, 1, 1], dim=0), dim=1).scalar()
    my_bottom = af.mean(af.mean(state.m[:, :, 0, 1], dim=0), dim=1).scalar()
    mx_a = mx_analytical(dz, Ms, extfield * Constants.mu0, RKKY_surface)
    my_a = my_analytical(dz, Ms, extfield * Constants.mu0, RKKY_surface)
    #mx_a_J = mx_analytical(dz, Js, extfield * Constants.mu0, RKKY_surface)
    #my_a_J = my_analytical(dz, Js, extfield * Constants.mu0, RKKY_surface)
    # TODO write only my_bottom
    #print(mx_a, my_a, mx_a_J, my_a_J)
    print(i, 'ext[T]={:2.3f}, mx={:1.3f}, my={:1.3f}, mz={:1.3f}, mxa={:1.3f}, mya={:1.3f}, mx_top={:1.3f}, mx_bottom={:1.3f}, my_top={:1.3f}, my_bottom={:1.3f}'.format(ext.h(state)[0, 0, 0, 1].scalar() * Constants.mu0, mx, my, mz, mx_a, my_a, mx_top, mx_bottom, my_top, my_bottom))
    stream.write("%e, %e, %e, %e, %e, %e, %e\n" %(extfield * Constants.mu0, mx, my, mz, mx_a, my_a, my_bottom))

stream.close()

# plotting data with gnuplot
f = open('plot.gpi', 'w')
f.write('set terminal pdf\n')
f.write('set output "' + sys.argv[1] + 'm.pdf"\n')
f.write('set xlabel "Hx [ns]"\n')
f.write('set ylabel "<m>"\n')
f.write('p "' + sys.argv[1] + '/m.dat" u 1:7 w l t "m_y"')
f.write(',"" u 1:6 w l t "m_{y, analytical}"')
f.close()

from os import system
system('gnuplot ' + sys.argv[1] +'plot.gpi')
system('evince '  + sys.argv[1] +'m.pdf')


## plotting data with gnuplot
#from os import system
#system('gnuplot -e "\
#    set terminal pdf;\
#    set output \'' + sys.argv[1] + 'm.pdf\';\
#    set xlabel \'Hx [ns]\';\
#    set ylabel \'<m>\';\
#    p \'' + sys.argv[1] + '/m.dat\' u 1:3 w l t \'<m_y>\',\
#    \'\' u 1:6 w l t \'<my_analytical>\',\
#    \'\' u 1:7 w l t \'my_bottom\';\
#"')
#    #\'\' u 1:6 w l t \'<my_analytical>\';\
#    #\'\' u 1:5 w l t \'<mx_analytical>\';\
#    #\'\' u 1:4 w l t \'<m_z>\';\
#    #\'\' u 1:3 w l t \'<m_y>\',\
#    #p \'' + sys.argv[1] + '/m.dat\' u 1:2 w l t \'<m_x>\',\
#
## show pdf with evince
#system('evince ' + sys.argv[1] +'m.pdf')

