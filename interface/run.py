import pth_mag
import ctypes
import arrayfire

arrayfire.info()
meshvar=pth_mag.pyMesh(3,4,5,2.,3.,4.)
print meshvar.n0()

#a = arrayfire.constant(2.,2,2,2,3)
#arrayfire.lock_array(a)
#print hex(ctypes.addressof(a.arr)),"hex(ctypes.addressof(a.arr))"
#meshvar.printme()
#x=test.pyTest()
#x.usearray(a)
