import arrayfire as af
import magnum_af
import sys
import os
import time

if sys.argv[1][-1] != "/":
    sys.argv[1] = sys.argv[1] + "/"
filepath = sys.argv[1]
#provided by bash #os.makedirs(filepath)
af.set_device(int(sys.argv[2]))
#af.set_backend("cpu")# TODO currently cpu backend segfaults when state object is created
af.info()


start = time.time()

# Physical dimensions in [m]
x = 5.e-7
y = 1.25e-7
z = 3.e-9
# Discretization 
nx = 100
ny = 25
nz = 1

# Creating objects
mesh=magnum_af.pyMesh(nx, ny, nz, x/nx, y/ny, z/nz)
m=af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)

param=magnum_af.pyParam()
#print (param.mu0)
param.ms = 8e5
param.A = 1.3e-11
param.alpha = 1.

m[1:-1,:,:,0] = af.constant(1.0, nx-2 ,ny, nz, 1, dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0, 1    ,ny, nz, 1, dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0, 1    ,ny, nz, 1, dtype=af.Dtype.f64);
state=magnum_af.pyState(mesh,param,m)
demag=magnum_af.pyDemagSolver(mesh,param)
exch=magnum_af.pyExchSolver(mesh,param)
Llg=magnum_af.pyLLG(demag,exch)

# Relaxing
print("relaxing 1ns")
stream = open(filepath+"m.dat", "w")
timer = time.time()
while state.t < 1e-9:
  Llg.llgstep(state)
  temp = state.m
  temp_mean = af.mean(af.mean(af.mean(temp, dim=0), dim=1), dim=2)
  #print(state.t(), temp_mean[0,0,0,0].scalar(), temp_mean[0,0,0,1].scalar(), temp_mean[0,0,0,2].scalar())
  stream.write("%e, %e, %e, %e\n" %(state.t, temp_mean[0,0,0,0].scalar(), temp_mean[0,0,0,1].scalar(), temp_mean[0,0,0,2].scalar()))
print("relaxed in", time.time() - timer, "[s]")

# Resetting alpha and adding Zeeman field
state.set_alpha(0.02)
zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/param.mu0
zeeswitch[0,0,0,1]=+4.3e-3/param.mu0
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch, nx, ny, nz)
zee=magnum_af.pyZee(zeeswitch)
Llg.add_terms(zee)

# Switching
print("switching 1ns")
timer = time.time()
while state.t < 2e-9:
  Llg.llgstep(state)
  temp = state.m
  temp_mean = af.mean(af.mean(af.mean(temp, dim=0), dim=1), dim=2)
  #print(state.t(), temp_mean[0,0,0,0].scalar(), temp_mean[0,0,0,1].scalar(), temp_mean[0,0,0,2].scalar())
  #print(state.meanxyz(0), state.meanxyz(1), state.meanxyz(2))
  stream.write("%e, %e, %e, %e\n" %(state.t, temp_mean[0,0,0,0].scalar(), temp_mean[0,0,0,1].scalar(), temp_mean[0,0,0,2].scalar()))
stream.close()
print("switched in", time.time() - timer, "[s]")
print("total time =", time.time() - start, "[s]")

#teststate=magnum_af.pyState(mesh,param,m) # testing wether teststate.m is correctly overwirtten with m_relax
#teststate.py_vti_reader("/home/paul/git/pth-mag/Data/Testing/py_interf/m_relax.vti")
#teststate.py_vti_writer_micro("/home/paul/git/pth-mag/Data/Testing/py_interf/m_reader")
