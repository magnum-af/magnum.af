import pth_mag
import arrayfire

arrayfire.info()
meshvar=pth_mag.pyMesh(15,15,15,1e-9,1e-9,1e-9)
print "n0= ", meshvar.n0()

param=pth_mag.pyParam()
param.print_gamma()
param.set_gamma(4e2)
param.set_ms(4e2)
param.print_gamma()
#TODO print param.D

#m=arrayfire.constant(42.0,100,100,100,3,dtype=arrayfire.Dtype.f64)
m=arrayfire.constant(1e19,15,15,15,3,dtype=arrayfire.Dtype.f64)
pystate=pth_mag.pyState(meshvar,param,m)
pystate.print_m()

demag=pth_mag.pyDemagSolver(meshvar,param)
demag.print_Nfft()
print "demag= ", demag.print_E(pystate)
print "demag_time= ",demag.cpu_time()

vec=pth_mag.pyVector(demag)
Llg=pth_mag.pyLLG(pystate,vec)
print "LLG= ", Llg.print_E(pystate)
print "LLG_cpu_time= ",Llg.cpu_time()

 
#TODO print pystate.get_m()

#a = arrayfire.constant(2.,2,2,2,3)
#arrayfire.lock_array(a)
#print hex(ctypes.addressof(a.arr)),"hex(ctypes.addressof(a.arr))"
#meshvar.printme()
#x=test.pyTest()
#x.usearray(a)
