import arrayfire as af
import interface as pth_mag
import ctypes
from copy import deepcopy

af.info()
meshvar=pth_mag.pyMesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)


param=pth_mag.pyParam()
param.gamma (2.211e5)
param.ms    (8e5)
param.A     (1.3e-11)
param.alpha (1)
#param.D_atom_axis_x(1)
#param.K_atom_axis(1.,0.,0.)
#param.K_atom_axis(100.,3.,30.)
#print "paramk= ",param.print_K_atom_axis_x() , param.print_K_atom_axis_y() , param.print_K_atom_axis_z()

m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
pystate=pth_mag.pyState(meshvar,param,m)
pystate.py_write_vtk()

##Alternative:
#m_test=af.Array()
#m_test_addr=pystate.get_m()
#m_test.arr = ctypes.c_void_p(m_test_addr)

#af.device.lock_array(m_test)


#af.device.lock_array(m_test)
#print "Test", af.mean(m_test)
#print af.device.is_locked_array(m_test)


demag=pth_mag.pyDemagSolver(meshvar,param)
exch=pth_mag.pyExchSolver(meshvar,param)
Llg=pth_mag.pyLLG(pystate,demag,exch)

print "relax --------------------"
print pystate.t()
#while pystate.t() < 1e-9:
Llg.llgstep(pystate)
print pystate.t()


print "switch --------------------"
Llg.set_state0_alpha(0.02)# this should be changed in cpp version

zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/param.print_mu0()
zeeswitch[0,0,0,1]=+4.3e-3/param.print_mu0()
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch,100,25,1)
zee=pth_mag.pyZee(zeeswitch,meshvar,param)
Llg.add_terms(zee)
#Llg2=pth_mag.pyLLG(pystate,demag,exch,zee)
print pystate.t()
#while pystate.t() < 2e-9:
Llg.llgstep(pystate)
print pystate.t()
print "end    --------------------"

#m_test=deepcopy(pystate.get_m())
#m_test=pystate.get_m()
print "af.mean(m_test)=", af.mean(pystate.get_m())
