
#HOTTIP
#https://stackoverflow.com/questions/44686590/initializing-cython-objects-with-existing-c-objects



#clib library_with_useful_functions
import ctypes
import arrayfire 
from libc.stdint cimport uintptr_t
from cython.operator cimport dereference as deref
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector

#from libc.stdint cimport uintptr_t

#from ctypes.wintypes import BOOL

#cdef extern from "<arrayfire.h>":
#  ctypedef void* af_array

#https://github.com/cython/cython/blob/master/Cython/Includes/libcpp/memory.pxd
#cdef extern from "<memory>" namespace "std":
#  cdef cppclass shared_ptr[T]:
#    shared_ptr()
#    shared_ptr(nullptr_t)
#    shared_ptr(T*)
#    shared_ptr(shared_ptr[T]&)
#    shared_ptr(shared_ptr[T]&, T*)
#    #shared_ptr(unique_ptr[T]&)
    
cdef extern from "../src/LLGTerm.hpp":
  cdef cppclass LLGTerm

cdef extern from "<arrayfire.h>" namespace "af":
  cdef cppclass array:
    array()

#cdef class pyVector:
#  pass; 
#cdef vector[Obj] vec;
#cdef  vec.push_back(Obj());

#cdef extern from "<vector>" namespace "std":
#    cdef cppclass vector[T]:
#        cppclass iterator:
#            T operator*()
#            iterator operator++()
#            bint operator==(iterator)
#            bint operator!=(iterator)
#        vector()
#        void push_back(T&)
#        T& operator[](int)
#        T& at(int)
#        iterator begin()
#        iterator end()
#cdef vector[int].iterator iter 

#cdef class pyVector:
#  cdef vector[shared_ptr]* thisptr
#  def __cinit__(self):
#    self.thisptr = new vector[shared_ptr] ()
#  def __dealloc__(self):
#    del self.thisptr

cdef class pyVector:
  cdef vector[shared_ptr[LLGTerm]]* thisptr
  def __cinit__(self):
    self.thisptr = new vector[shared_ptr[LLGTerm]] ()
  def __dealloc__(self):
    del self.thisptr
  #def push_back(self, shared_ptr[LLGTerm] a):
  #  self.thisptr.push_back(a)

cdef extern from "../src/llg.hpp":
  cdef cppclass LLG:
    LLG (State state_in, vector[shared_ptr[LLGTerm]] vector_in);
    array llgstep(State& state);
    double E(const State& state);
    double cpu_time();

cdef class pyLLG:
  cdef LLG* thisptr
  def __cinit__(self, pyState state_in, pyVector vector_in):
    self.thisptr = new LLG (deref(state_in.thisptr), deref(vector_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
#TODO  def llgstep(self, pyState state_in):
#    return self.thisptr.llgstep(deref(state_in.thisptr))
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.cpu_time()

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
  def set_ms(self,ms):
    self.thisptr.ms=ms
  
  def print_gamma(self):
    print self.thisptr.gamma


cdef extern from "../src/state.hpp":
  cdef cppclass State:
    State (Mesh mesh_in, Param param_in, long int m_in);
    void print_m();
    long int get_m();

cdef class pyState:
  cdef State* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in, m_in):
    self.thisptr = new State (deref(mesh_in.thisptr), deref(param_in.thisptr), ctypes.addressof(m_in.arr))  
  def __dealloc__(self):
    del self.thisptr
  def print_m(self):
    self.thisptr.print_m()

#  #TODO
#  def get_m(self):
#    cdef void **a = (void **) self.thisptr.get_m()
#    cdef array m_out = * ( array (* a))
#    return m_out
#    #return self.thisptr.get_m()


cdef extern from "../src/micro_demag.hpp":
  cdef cppclass DemagSolver:
    DemagSolver (Mesh mesh_in, Param param_in);
    double E(const State& state);
    void print_Nfft();
    double get_cpu_time();

cdef class pyDemagSolver:
  cdef DemagSolver* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    self.thisptr = new DemagSolver (deref(mesh_in.thisptr), deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
  def print_Nfft(self):
    self.thisptr.print_Nfft()
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()






#cdef uintptr_t adr = <uintptr_t>ctypes.addressof(foo_ct.contents)
#cy_use_struct(<Foo*>adr)

#cdef extern from "teststate.hpp":
#  cdef cppclass testState:
#    testState (Mesh mesh_in, Param param_in, long int a);
#    void printn0();
#
#cdef class pytestState:
#  cdef testState* thisptr
#  def printn0(self):
#    return self.thisptr.printn0()
# 
#  def __cinit__(self, pyMesh mesh_in, pyParam param_in, a_in):
#    self.thisptr = new testState (deref(mesh_in.thisptr), deref(param_in.thisptr),ctypes.addressof(a_in.arr))  
#  #def __cinit__(self, pyMesh mesh_in, pyParam param_in):
#  #  self.thisptr = new testState (deref(mesh_in.thisptr), deref(param_in.thisptr))  
#  def __dealloc__(self):
#    del self.thisptr


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

#cdef extern from "teststate.hpp":
#  cdef cppclass testState:
#    testState (Mesh mesh_in, Param param_in);
#    void printn0();
#
##cdef uintptr_t adr = <uintptr_t>ctypes.addressof(foo_ct.contents)
##cy_use_struct(<Foo*>adr)
#
#cdef class pytestState:
#  cdef testState* thisptr
#  #cdef pyMesh mesh
#  #cdef pyParam param
#  def printn0(self):
#    return self.thisptr.printn0()
# 
#  #def __cinit__(self, mesh_in,  param_in):
#  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
#    #self.thisptr = new testState (<uintptr_t>ctypes.addressof(mesh_in),<uintptr_t>ctypes.addressof(param_in))  
#    mesh=mesh_in
#    param=param_in
#    #cdef uintptr_t adr1 = <uintptr_t>ctypes.addressof(mesh_in)
#    #cdef uintptr_t adr2 = <uintptr_t>ctypes.addressof(param_in)
#    #self.thisptr = new testState (<Mesh*>adr1,<Param*>adr2)  
#    self.thisptr = new testState (deref(mesh_in.thisptr), deref(param_in.thisptr))  
#  def __dealloc__(self):
#    del self.thisptr
#
