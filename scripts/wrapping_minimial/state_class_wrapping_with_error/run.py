import ctypes
import arrayfire
import test

x = test.pyTest()
a = arrayfire.constant(2., 4, dtype=arrayfire.Dtype.f64)
print "py: a=", a

x.py_initialize_m(a)
x.py_manipulate_m()

b = x.py_get_m()
print "b =", b

x.py_manipulate_m()
c = x.py_get_m()
print "c =", c
