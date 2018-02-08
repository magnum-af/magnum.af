import arrayfire as af
import pth_mag
import ctypes
from copy import deepcopy

#af.set_backend("cpu")
af.info()
meshvar=pth_mag.pyMesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)


param=pth_mag.pyParam()
param.ms    (8e5)
param.A     (1.3e-11)
param.alpha (1)

m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
pystate=pth_mag.pyState(meshvar,param,m)

demag=pth_mag.pyDemagSolver(meshvar,param)
exch=pth_mag.pyExchSolver(meshvar,param)
Llg=pth_mag.pyLLG(pystate,demag,exch)

print "relax --------------------"
print pystate.t()
while pystate.t() < 1e-9:
  Llg.llgstep(pystate)
print pystate.t()
pystate.py_vti_writer_micro("/home/pth/git/pth-mag/Data/Testing/py_interf/m_test")

print "switch --------------------"
Llg.set_state0_alpha(0.02)# this should be changed in cpp version

zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/param.print_mu0()
zeeswitch[0,0,0,1]=+4.3e-3/param.print_mu0()
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch,100,25,1)
zee=pth_mag.pyZee(zeeswitch,meshvar,param)
Llg.add_terms(zee)
print pystate.t()
while pystate.t() < 2e-9:
  Llg.llgstep(pystate)
print pystate.t()
print "end    --------------------"
print "af.mean(m_test)=", af.mean(pystate.get_m())
