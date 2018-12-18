import ctypes
import arrayfire 

cdef extern from "source.hpp":
  cdef cppclass Wrap:
    void c_py_to_cpp(long int a_addr)
    void c_calc_B()
    long int c_cpp_to_py()

cdef class pyWrap:
  cdef Wrap* thisptr
  def __cinit__(self):
    self.thisptr = new Wrap()
  def __dealloc__(self):
    del self.thisptr
  def py_to_cpp(self, a):
    self.thisptr.c_py_to_cpp(ctypes.addressof(a.arr))
  def cpp_to_py(self):
    return self.thisptr.c_cpp_to_py()
  def calc_B(self):
    self.thisptr.c_calc_B()
