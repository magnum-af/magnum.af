import arrayfire as af
import magnum_af
import ctypes
from copy import deepcopy

#af.set_backend("cpu")
af.info()
meshvar=magnum_af.pyMesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)


param=magnum_af.pyParam()
param.ms    (8e5)
param.A     (1.3e-11)
param.alpha (1)

m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
pystate=magnum_af.pyState(meshvar,param,m)

demag=magnum_af.pyDemagSolver(meshvar,param)
exch=magnum_af.pyExchSolver(meshvar,param)
Llg=magnum_af.pyLLG(pystate,demag,exch)

print "relax --------------------"
print pystate.t()
while pystate.t() < 1e-9:
  time1=pystate.t()
  Llg.llgstep(pystate)
  time2=pystate.t()
  print Llg.print_stepsize(), time2-time1, time1, time2, pystate.t()
print pystate.t()
#pystate.py_vti_writer_micro("/home/paul/git/pth-mag/Data/Testing/py_interf/m_relax")

#teststate=magnum_af.pyState(meshvar,param,m) # testing wether teststate.m is correctly overwirtten with m_relax
#teststate.py_vti_reader("/home/paul/git/pth-mag/Data/Testing/py_interf/m_relax.vti")
#teststate.py_vti_writer_micro("/home/paul/git/pth-mag/Data/Testing/py_interf/m_reader")

print "switch --------------------"
Llg.set_state0_alpha(0.02)# this should be changed in cpp version

zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/param.print_mu0()
zeeswitch[0,0,0,1]=+4.3e-3/param.print_mu0()
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch,100,25,1)
zee=magnum_af.pyZee(zeeswitch)
Llg.add_terms(zee)
print pystate.t()
intx=0
inty=0
intz=0
while pystate.t() < 2e-9:
  time1=pystate.t()
  Llg.llgstep(pystate)
  time2=pystate.t()
  intx+=pystate.meanxyz(0)*Llg.print_stepsize()
  inty+=pystate.meanxyz(1)*Llg.print_stepsize()
  intz+=pystate.meanxyz(2)*Llg.print_stepsize()
  #print Llg.print_stepsize(), time2-time1, time1, time2, pystate.t()
  #print pystate.meanxyz(0), pystate.meanxyz(1), pystate.meanxyz(2)
  #  std::cout << "cpp h= "<< h << "t2-t1 "<< time2-time1 << " time1 "<< time1<< " time2 " << time2 << "state.t "<< state.t << std::endl;
  #  std::cout.precision(12);
  print intx, inty, intz
print pystate.t()
print "finished  --------------------"
#TODO second call causes segfault:
print "af.mean(m_test)=", af.mean(pystate.get_m())
print "af.mean(m_test)=", af.mean(pystate.get_m())
