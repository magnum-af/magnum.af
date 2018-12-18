import ctypes
import arrayfire 
from libc.stdint cimport uintptr_t

cdef extern from "<arrayfire.h>":
  ctypedef void* af_array

cdef extern from "<arrayfire.h>" namespace "af":
  cdef cppclass array:
    array()

cdef extern from "cpp_test.h":
  cdef cppclass Test:
    Test();
    void init_m(long int aptr);
    void print_m();
    long int get_m();

cdef class pyTest:
  cdef Test* thisptr # hold a C++ instance
  def __cinit__(self):
    self.thisptr = new Test()

  # Leaving this uncommented leads to error in cleanup
  # Commenting this out allows py_get_m() being called once without error, second call yields double-free or corruption
  #def __dealloc__(self):
  #  del self.thisptr

  def py_init_m(self, a):
    self.thisptr.init_m(ctypes.addressof(a.arr))

  def py_print_m(self):
    self.thisptr.print_m()

  def py_get_m(self):
    b_adr = self.thisptr.get_m()    
    b = arrayfire.Array()
    b.arr = ctypes.c_void_p(b_adr)
    return b
