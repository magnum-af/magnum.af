import pth_mag
import arrayfire as af

af.info()
meshvar=pth_mag.pyMesh(  100,25,1,5.e-7/100,1.25e-7/25,3.e-9)
m=af.constant(0.0,100,25,1,3,dtype=af.Dtype.f64)
print "n0= ", meshvar.n0()

param=pth_mag.pyParam()
param.gamma (2.211e5)
param.ms    (8e5)
param.A     (1.3e-11)
param.alpha (1)
#TODO param.D_axis(3)
#TODO print param.D

#m=arrayfire.constant(42.0,100,100,100,3,dtype=arrayfire.Dtype.f64)
m[1:-1,:,:,0] = af.constant(1.0,100-2,25,1,1,dtype=af.Dtype.f64);
m[0,:,:,1]    = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
m[-1,:,:,1]   = af.constant(1.0,1    ,25,1,1,dtype=af.Dtype.f64);
pystate=pth_mag.pyState(meshvar,param,m)
#pystate.print_m()

demag=pth_mag.pyDemagSolver(meshvar,param)
exch=pth_mag.pyExchSolver(meshvar,param)
#demag2=pth_mag.pyDemagSolver(meshvar,param)
#demag3=pth_mag.pyDemagSolver(meshvar,param)
#demag.print_Nfft()
#print "demag= ", demag.print_E(pystate)
#print "demag_time= ",demag.cpu_time()

#vec=pth_mag.pyVector(demag)
#Llg=pth_mag.pyLLG(pystate,vec)
#print demag.addr()
Llg=pth_mag.pyLLG(pystate,demag,exch)
#Llg=pth_mag.pyLLG(pystate,demag,demag2,demag3)
#print "LLG= ", Llg.print_E(pystate)
#print "LLG_cpu_time= ",Llg.cpu_time()
print "relax --------------------"
print pystate.t()
#while pystate.t() < 1e-9:
Llg.llgstep(pystate)
print pystate.t()


print "switch --------------------"
Llg.set_state0_alpha(0.02)# this should be changed in cpp version

print "param =", param.mu0()
zeeswitch = af.constant(0.0,1,1,1,3,dtype=af.Dtype.f64)
zeeswitch[0,0,0,0]=-24.6e-3/param.mu0()
zeeswitch[0,0,0,1]=+4.3e-3/param.mu0()
zeeswitch[0,0,0,2]=0.0
zeeswitch = af.tile(zeeswitch,100,25,1)
zee=pth_mag.pyZee(zeeswitch,meshvar,param)
Llg2=pth_mag.pyLLG(pystate,demag,exch,zee)
print pystate.t()
#while pystate.t() < 2e-9:
Llg2.llgstep(pystate)
print pystate.t()
 
#TODO print pystate.get_m()

#a = arrayfire.constant(2.,2,2,2,3)
#arrayfire.lock_array(a)
#print hex(ctypes.addressof(a.arr)),"hex(ctypes.addressof(a.arr))"
#meshvar.printme()
#x=test.pyTest()
#x.usearray(a)
