import ctypes
import wrap
import arrayfire as af

wrapper=wrap.pyWrap()

A=af.constant(1.0,2,3,2,3, af.Dtype.f64)

wrapper.py_to_cpp(A)
wrapper.calc_B()
B_addr = wrapper.cpp_to_py()


B = af.Array(B_addr,A.shape, af.Dtype.f64) #This line produces wrong values and needs fixing
print "py: B = ", B
