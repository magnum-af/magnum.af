import ctypes
import arrayfire
import test

x = test.pyTest()
a = arrayfire.randu(4)
print "py: a=", a

x.py_init_m(a)
x.py_print_m()
b = x.py_get_m()
print "b =", b
#c = x.py_get_m()
#print "c =", c
