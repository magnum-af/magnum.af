
#HOTTIP
#https://stackoverflow.com/questions/44396749/handling-c-arrays-in-cython-with-numpy-and-pytorch

#TODO have a look on
#from libcpp.memory cimport shared_ptr
#
#cdef class Holder:
#    cdef shared_ptr[cpp_class] ptr
#
#    @staticmethod
#    cdef make_holder(shared_ptr[cpp_class] ptr):
#       cdef holder = Holder() # empty class
#       holder.ptr = ptr
#       return holder

#https://stackoverflow.com/questions/44686590/initializing-cython-objects-with-existing-c-objects
#https://stackoverflow.com/questions/47044866/how-to-convert-python-object-to-c-type-in-cython



#clib library_with_useful_functions
import ctypes
import arrayfire 
from libc.stdint cimport uintptr_t
from cython.operator cimport dereference as deref
from libcpp.memory cimport shared_ptr#, make_shared
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
#  cdef vector* thisptr
#  def __cinit__(self):
#    self.thisptr = new vector ()
#  def __dealloc__(self):
#    del self.thisptr

#cdef class pyVector:
#  cdef vector[shared_ptr[LLGTerm]]* thisptr
#  def __cinit__(self, shared_ptr[LLGTerm] T_in):
#    self.thisptr = new vector[shared_ptr[LLGTerm]] (T_in)
#  def __dealloc__(self):
#    del self.thisptr
##  def py_push_back(self,  T):

cdef class pyVector:
  cdef vector[shared_ptr[LLGTerm]]* thisptr
  def __cinit__(self):
    self.thisptr = new vector[shared_ptr[LLGTerm]] ()
  def __dealloc__(self):
    del self.thisptr
#  def push_back(self,  shared_ptr[LLGTerm] a):
#    self.thisptr.push_back( a)

cdef extern from "<utility>" namespace "std" nogil:
    #cdef shared_ptr[LLGTerm] move(unique_ptr[LLGTerm])
    cdef shared_ptr[LLGTerm] move(shared_ptr[LLGTerm])


#cdef public newDoodadFromSptr(shared_ptr[LLGTerm] _d):
#    d = LLGTerm(init=False)
#    d.thisptr = move(_d)
#
#    return d

#https://dmtn-013.lsst.io/
#cpdef as_list(self):
#    cdef vector[shared_ptr[_Doodad]] v = self.inst.as_vector()
#
#    results = []
#    for item in v:
#        d = Doodad(init=False)
#        d.thisptr = move(item)
#        results.append(d)
#
#    return results

cdef extern from "../src/llg.hpp":
  cdef cppclass LLG:
    LLG (State state_in, vector[shared_ptr[LLGTerm]] vector_in);
    array llgstep(State& state);
    double E(const State& state);
    double cpu_time();

cdef class pyLLG:
  cdef LLG* thisptr
  #def __cinit__(self, pyState state_in, pyVector vector_in):
  def __cinit__(self, pyState state_in, *args):
    cdef vector[shared_ptr[LLGTerm]] vector_in
    #cdef shared_ptr[LLGTerm] temp = shared_ptr[LLGTerm] (<LLGTerm*><size_t>terms.pythisptr())
    #cdef shared_ptr[LLGTerm] temp = shared_ptr[LLGTerm] (<LLGTerm*><size_t>terms.addr())
    #cdef shared_ptr[LLGTerm] temp = shared_ptr[LLGTerm] (<LLGTerm*><size_t>ctypes.addressof(terms))
    #vector_in.push_back(temp)
    for arg in args:
      vector_in.push_back(shared_ptr[LLGTerm] (<LLGTerm*><size_t>arg.pythisptr()))
    #for term in terms:
    #vector_in.push_back(shared_ptr[LLGTerm](<LLGTerm*>terms.thisptr))   
    #vector_in.push_back[shared_ptr[LLGTerm]](terms)
    #cdef shared_ptr[LLGTerm] cterm = make_shared[LLGTerm](terms.thisptr)
    #vector_in.push_back(make_shared[LLGTerm](deref(term.thisptr)))
    self.thisptr = new LLG (deref(state_in.thisptr), vector_in)  
    #self.thisptr = new LLG (deref(state_in.thisptr), <vector[shared_ptr[LLGTerm]]>vector_in))  
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

cdef class pyMesh:
  cdef Mesh* thisptr
  def __cinit__(self,int a, int b, int c, double d, double e, double f):
    self.thisptr = new Mesh(a,b,c,d,e,f)

  def __dealloc__(self):
    del self.thisptr

  def n0(self):
    return self.thisptr.n0

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

cdef class pyParam:
  cdef Param* thisptr
  def __cinit__(self):
    self.thisptr = new Param ()  
  def __dealloc__(self):
    del self.thisptr

  def gamma(self,value):
    self.thisptr.gamma=value
  def ms(self,value):
    self.thisptr.ms=value
  def alpha(self,value):
    self.thisptr.alpha=value
  def A(self,value):
    self.thisptr.A=value
  def p(self,value):
    self.thisptr.p=value
  def mode(self,value):
    self.thisptr.mode=value
  def D(self,value):
    self.thisptr.D=value
  def Ku1(self,value):
    self.thisptr.Ku1=value
  def J_atom(self,value):
    self.thisptr.J_atom=value
  def D_atom(self,value):
    self.thisptr.D_atom=value
  def K_atom(self,value):
    self.thisptr.K_atom=value
  #TODO
  #uudef D_axis(self,value):
  #  self.thisptr.D_axis[7]=value
  
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
  #TODO this causes double free coruption!!!!
  #def __dealloc__(self):
  #  del self.thisptr
  def print_Nfft(self):
    self.thisptr.print_Nfft()
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr






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
