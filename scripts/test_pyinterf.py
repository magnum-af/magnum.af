import arrayfire as af
import magnumaf
import ctypes
from copy import deepcopy

#af.set_backend("cpu")
af.info()
meshvar=magnumaf.Mesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)


material=magnumaf.Material()
state.Ms    (8e5)
material.A     (1.3e-11)
material.alpha (1)

m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
pystate=magnumaf.State(meshvar,material,m)

demag=magnumaf.DemagField(meshvar,material)
exch=magnumaf.ExchangeField(meshvar,material)
Llg=magnumaf.LLGIntegrator([pystate,demag,exch])

print "relax --------------------"
print pystate.t()
while pystate.t() < 1e-9:
  Llg.step(pystate)
print pystate.t()
pystate.write_vti("/home/pth/git/magnum.af/Data/Testing/py_interf/m_relax")

teststate=magnumaf.State(meshvar,material,m) # testing wether teststate.m is correctly overwirtten with m_relax
teststate.read_vti("/home/pth/git/magnum.af/Data/Testing/py_interf/m_relax.vti")
teststate.write_vti("/home/pth/git/magnum.af/Data/Testing/py_interf/m_reader")

print "switch --------------------"
Llg.set_state0_alpha(0.02)# this should be changed in cpp version

zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/material.print_mu0()
zeeswitch[0,0,0,1]=+4.3e-3/material.print_mu0()
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch,100,25,1)
zee=magnumaf.ExternalField(zeeswitch)
Llg.add_terms(zee)
print pystate.t()
while pystate.t() < 2e-9:
  Llg.step(pystate)
print pystate.t()
print "end    --------------------"
print "af.mean(m_test)=", af.mean(pystate.get_m())
