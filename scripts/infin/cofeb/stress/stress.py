import arrayfire as af
import numpy as np
from magnum_af import *
import time

def ellipse(n0, n1, n2):
    m = af.constant(0.0, n0, n1, n2, 3, dtype=af.Dtype.f64)
    for ix in range (0, n0):
        for iy in range(0, n1):
            a= n0/2
            b= n1/2
            rx=ix-n0/2.
            ry=iy-n1/2.
            r = pow(rx,2)/pow(a,2)+pow(ry,2)/pow(b,2)
            if(r<1):
                m[ix,iy,:,0]=1
    return m

def np_ellipse(n0, n1, n2):
    m = np.zeros((n0, n1, n2, 3));
    for ix in range (0, n0):
        for iy in range(0, n1):
            a= n0/2
            b= n1/2
            rx=ix-n0/2.
            ry=iy-n1/2.
            r = pow(rx,2)/pow(a,2)+pow(ry,2)/pow(b,2);
            if(r<1):
                m[ix,iy,:,0]=1
    return af.from_ndarray(m)

af.info()
filepath = "todel/pyrun1/"
x = 800e-9
y = 800e-9
z = 1.3e-3/1.056e6
nx = 250
ny = 250
nz = 1
mesh=pyMesh(nx, ny, nz, x/nx, y/ny, z/nz)
m=af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)

#print (ellipse(nx, ny, nz))

param=pyParam()
param.ms(1.58/param.print_mu0())
param.A(15e-12)
param.Ku1(1.3e-3/z)

param_stress=pyParam()
param_stress.ms(1.58/param.print_mu0())
param_stress.A(15e-12)
param_stress.Ku1(0.1* 1.3e-3/z) # TODO correct value
param_stress.Ku1_axis(1, 0, 0) 
print ("Ku1 axis =", param_stress.print_Ku1_axis())

start = time.time()
pystate=pyState(mesh, param, np_ellipse(nx, ny, nz))
print ("Init ellipse [s]= ", time.time() - start)

start = time.time()
demag=pyDemagSolver(mesh, param)
exch=pyExchSolver(mesh, param)
aniso_z = pyMicroAniso(mesh, param)
aniso_stress = pyMicroAniso(mesh, param_stress)

zeeswitch = af.constant(0.0, nx, ny, nz, 3,dtype=af.Dtype.f64)
zee = pyZee(zeeswitch)
zee.set_xyz(pystate, 50e-3/param.print_mu0(), 0, 0)
print ("Init terms [s]= ", time.time() - start)

#minimizer = pyLbfgsMinimizer(demag, exch, aniso_z, aniso_stress, zee)
#we want somthing like this: minimizer.zee.set_xyz(pystate, 50e-3/param.print_mu0(), 0, 0)

#print ("GetTimeCalcHeff() [s]= ", minimizer.pyGetTimeCalcHeff())

#pystate.py_vti_writer_micro("todel/py_test_m_init")
#start = time.time()
#minimizer.pyMinimize(pystate)
#print ("Minimize [s]= ", time.time() - start)
#pystate.py_vti_writer_micro("todel/py_test_m_minimized_0")


A = 0.05/param.print_mu0()
steps = 100
print ("A= ", A)
for i in range(0, steps):
    phi = 2. * np.pi * i/steps;
    print ("phi= ", phi)
    zee.set_xyz(pystate, A * np.cos(phi), A * np.sin(phi), 0)
    minimizer = pyLbfgsMinimizer(demag, exch, aniso_z, aniso_stress, zee)
    start = time.time()
    minimizer.pyMinimize(pystate)
    print ("Minimize " + str(i) + " in [s]= ", time.time() - start)
    pystate.py_vti_writer_micro(filepath + "m_"+ str(i))

#minimizer=[]
#A = 0.05/param.print_mu0()
#steps = 100
#print ("A= ", A)
#for i in range(0, steps):
#    phi = 2. * np.pi * i/steps;
#    print ("phi= ", phi)
#    zee.set_xyz(pystate, A * np.cos(phi), A * np.sin(phi), 0)
#    minimizer.append (pyLbfgsMinimizer(demag, exch, aniso_z, aniso_stress, zee))
#    start = time.time()
#    minimizer[i].pyMinimize(pystate)
#    print ("Minimize " + str(i) + " in [s]= ", time.time() - start)
#    pystate.py_vti_writer_micro(filepath + "m_"+ str(i))

################################################
# TODO this is causing segfault:
#minimizer.delete_last_term()# TODO this causes segfault!!!!!!
#zee.set_xyz(pystate, 0, 50e-3/param.print_mu0(), 0)
#minimizer.add_terms(zee)
#start = time.time()
#print ("test0")
#minimizer.pyMinimize(pystate)
#print ("Minimize [s]= ", time.time() - start)
