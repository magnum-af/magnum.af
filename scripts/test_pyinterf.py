import arrayfire as af
import magnum_af
import ctypes
from copy import deepcopy

#af.set_backend("cpu")
af.info()
meshvar=magnum_af.Mesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)


param=magnum_af.Param()
param.ms    (8e5)
param.A     (1.3e-11)
param.alpha (1)

m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
pystate=magnum_af.State(meshvar,param,m)

demag=magnum_af.DemagSolver(meshvar,param)
exch=magnum_af.ExchSolver(meshvar,param)
Llg=magnum_af.NewLlg(pystate,demag,exch)

print "relax --------------------"
print pystate.t()
while pystate.t() < 1e-9:
  Llg.llgstep(pystate)
print pystate.t()
pystate.py_vti_writer_micro("/home/pth/git/magnum.af/Data/Testing/py_interf/m_relax")

teststate=magnum_af.State(meshvar,param,m) # testing wether teststate.m is correctly overwirtten with m_relax
teststate.py_vti_reader("/home/pth/git/magnum.af/Data/Testing/py_interf/m_relax.vti")
teststate.py_vti_writer_micro("/home/pth/git/magnum.af/Data/Testing/py_interf/m_reader")

print "switch --------------------"
Llg.set_state0_alpha(0.02)# this should be changed in cpp version

zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/param.print_mu0()
zeeswitch[0,0,0,1]=+4.3e-3/param.print_mu0()
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch,100,25,1)
zee=magnum_af.Zee(zeeswitch)
Llg.add_terms(zee)
print pystate.t()
while pystate.t() < 2e-9:
  Llg.llgstep(pystate)
print pystate.t()
print "end    --------------------"
print "af.mean(m_test)=", af.mean(pystate.get_m())
