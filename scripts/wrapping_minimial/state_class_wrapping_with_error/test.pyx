import arrayfire 
import ctypes

cdef extern from "<arrayfire.h>":
  ctypedef void* af_array

cdef extern from "<arrayfire.h>" namespace "af":
  cdef cppclass array:
    array()

cdef extern from "cpp_test.h":
  cdef cppclass Test:
    Test();
    void initialize_m(long int aptr);
    void manipulate_m();
    long int get_m();

cdef class pyTest:
  cdef Test* thisptr # hold a C++ instance
  def __cinit__(self):
    self.thisptr = new Test()

  def __dealloc__(self):
    del self.thisptr

  def py_initialize_m(self, a):
    self.thisptr.initialize_m(ctypes.addressof(a.arr))

  def py_manipulate_m(self):
    self.thisptr.manipulate_m()

  def py_get_m(self):
    b_adr = self.thisptr.get_m()    
    b = arrayfire.Array()
    b.arr = ctypes.c_void_p(b_adr)
    return b
