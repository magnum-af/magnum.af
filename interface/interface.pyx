import ctypes
import arrayfire 
from libc.stdint cimport uintptr_t
from ctypes.wintypes import BOOL

cdef extern from "<arrayfire.h>":
  ctypedef void* af_array

cdef extern from "<arrayfire.h>" namespace "af":
  cdef cppclass array:
    array()

cdef extern from "../src/mesh.hpp":
  cdef cppclass Mesh:
    int n0,n1,n2;
    double dx,dy,dz;
    int n0_exp, n1_exp, n2_exp;
    Mesh (int, int, int, double, double, double)
    void printme()

cdef extern from "../src/param.hpp":
  cdef cppclass Param:
#    double mu0;
    double gamma,ms,A,alpha;
    double p;
#    bool afsync;
    int mode;
    double D;
    double Ku1;
    double D_axis[3];
    double Ku1_axis[3];
    double J_atom;
    double D_atom;
    double D_atom_axis[3];
    double K_atom;
    double K_atom_axis[3];
 

#cdef class pyTest:
#  cdef Test* thisptr # hold a C++ instance
#  def __cinit__(self):
#    self.thisptr = new Test()
#
#  def __dealloc__(self):
#    del self.thisptr
#
#  def usearray(self, a):
#    adr = ctypes.addressof(a.arr)
#    print hex(adr), "hex(adr) in usearray";
#    arrayfire.info()
#    self.thisptr.usearray(adr)

cdef class pyMesh:
  cdef Mesh* thisptr
  def __cinit__(self,int a, int b, int c, double d, double e, double f):
    self.thisptr = new Mesh(a,b,c,d,e,f)

  def __dealloc__(self):
    del self.thisptr

  def n0(self):
    return self.thisptr.n0


#  def printme(self):
#    self.thisptr.printme()
    
