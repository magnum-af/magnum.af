import pth_mag
import ctypes
import arrayfire

arrayfire.info()
meshvar=pth_mag.pyMesh(42,4,5,2.,3.,4.)
print "n0= ", meshvar.n0()

param=pth_mag.pyParam()
param.print_gamma()
param.set_gamma(4e2)
param.print_gamma()
#TODO print param.D

m=arrayfire.constant(42.0,3,3,3,1,dtype=arrayfire.Dtype.f64)
pystate=pth_mag.pytestState(meshvar,param,m)
pystate.printn0()

#a = arrayfire.constant(2.,2,2,2,3)
#arrayfire.lock_array(a)
#print hex(ctypes.addressof(a.arr)),"hex(ctypes.addressof(a.arr))"
#meshvar.printme()
#x=test.pyTest()
#x.usearray(a)
