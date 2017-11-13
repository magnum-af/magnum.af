import ctypes
import arrayfire 
from libc.stdint cimport uintptr_t
from cython.operator cimport dereference as deref
from libc.stdint cimport uintptr_t

#from ctypes.wintypes import BOOL

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

cdef class pyMesh:
  cdef Mesh* thisptr
  def __cinit__(self,int a, int b, int c, double d, double e, double f):
    self.thisptr = new Mesh(a,b,c,d,e,f)

  def __dealloc__(self):
    del self.thisptr

  def n0(self):
    return self.thisptr.n0

cdef class pyParam:
  cdef Param* thisptr
  def __cinit__(self):
    self.thisptr = new Param ()  

  def __dealloc__(self):
    del self.thisptr

  def set_gamma(self,gamma_in):
    self.thisptr.gamma=gamma_in
  
  def print_gamma(self):
    print self.thisptr.gamma


##TODO
#cdef extern from "../src/state.hpp":
#  cdef cppclass State:
#    State (Mesh mesh_in, Param param_in, array m_in);
#
#cdef class pyState:
#  cdef State* thisptr
#  #cdef pyMesh mesh
#  #cdef pyParam param
#  def __cinit__(self, pyMesh mesh_in, pyParam param_in, array m_in):
#    self.thisptr = new State (deref(mesh_in.thisptr), deref(param_in.thisptr), ctypes.addressof(m_in.arr))  
#  def __dealloc__(self):
#    del self.c_state
##TODO END

cdef extern from "teststate.hpp":
  cdef cppclass testState:
    testState (Mesh mesh_in, Param param_in);
    void printn0();

#cdef uintptr_t adr = <uintptr_t>ctypes.addressof(foo_ct.contents)
#cy_use_struct(<Foo*>adr)

cdef class pyState:
  cdef testState* thisptr
  cdef pyMesh mesh
  cdef pyParam param
  def printn0(self):
    return self.thisptr.printn0()
 
  #def __cinit__(self, mesh_in,  param_in):
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    #self.thisptr = new testState (<uintptr_t>ctypes.addressof(mesh_in),<uintptr_t>ctypes.addressof(param_in))  
    mesh=mesh_in
    param=param_in
    #cdef uintptr_t adr1 = <uintptr_t>ctypes.addressof(mesh_in)
    #cdef uintptr_t adr2 = <uintptr_t>ctypes.addressof(param_in)
    #self.thisptr = new testState (<Mesh*>adr1,<Param*>adr2)  
    self.thisptr = new testState (deref(mesh_in.thisptr), deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr


#cdef class pyState:
#  cdef State* thisptr
#  cdef pyMesh mesh
#  cdef pyState state
#  def __cinit__(self, pyMesh mesh_in, pyParam param_in, array m_in):
#    self.thisptr = new State (ctypes.addressof(mesh_in.thisptr), param_in, m_in)  
#
#  def __dealloc__(self):
#    del self.thisptr

  
  #def __cinit__(self, Mesh mesh_in, Param param_in, array m_in):
  #def __cinit__(self, pyMesh mesh_in, pyParam param_in, array m_in):

  #def __cinit__(self, mesh_in, param_in, m_in):
  #  self.thisptr = new State (mesh_in, param_in, m_in)  

  #def __dealloc__(self):
  #  del self.thisptr

  
#  def printme(self):
#    self.thisptr.printme()
    

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

