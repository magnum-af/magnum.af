#!python
#distutils: language = c++
#cython: language_level=3

#HOTTIPS
#https://stackoverflow.com/questions/44396749/handling-c-arrays-in-cython-with-numpy-and-pytorch
#https://stackoverflow.com/questions/44686590/initializing-cython-objects-with-existing-c-objects
#https://stackoverflow.com/questions/47044866/how-to-convert-python-object-to-c-type-in-cython
#clib library_with_useful_functions

## numpy to arrayfire
#af.interop.np_to_af_array
## arrayfire to numpy
#af.Array.__array__()

import arrayfire as af 
import copy
from libc.stdint cimport uintptr_t
from libcpp.memory cimport shared_ptr#, make_shared
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref
import ctypes
from math import sqrt


cdef extern from "../../src/llg_terms/LLGTerm.hpp":
  cdef cppclass LLGTerm

cdef extern from "<arrayfire.h>" namespace "af":
  cdef cppclass array:
    array()

cdef extern from "../../src/mesh.hpp":
#cdef extern from "/home/pth/git/cmakecython/src/mesh.hpp":
  cdef cppclass Mesh:
    int n0,n1,n2;
    double dx,dy,dz;
    int n0_exp, n1_exp, n2_exp;
    Mesh (int, int, int, double, double, double)


#NOTE#@cython.embedsignature(True)# error: Cdef functions/classes cannot take arbitrary decorators. https://stackoverflow.com/questions/42668252/cython-cdef-class-not-displaying-doc-string-or-init-parameters
# Docstring does work, todo: check type etc. 
cdef class pyMesh:
  """
  Class defining number of nodes and discretization.
  """
  def __init__(self, n0, n1, n2, dx, dy, dz): # todo: For dockstring, but currently does not change signature
    pass

  cdef Mesh* thisptr
  def __cinit__(self,int n0, int n1, int n2, double dx, double dy, double dz):
    self.thisptr = new Mesh(n0, n1, n2, dx, dy, dz)
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def n0(self):
    return self.thisptr.n0
  def print_n0(self):
    print (self.thisptr.n0)


cdef extern from "../../src/state.hpp":
  cdef cppclass State:
    State (Mesh mesh_in, Param param_in, long int m_in);
    State (Mesh mesh_in, Param param_in, long int aptr, long int evaluate_mean_ptr);
    void set_m(long int aptr);
    double t;
    array m;
    Mesh mesh;
    Param param;
    long int get_m_addr();

    void _vti_writer_micro(string outputname);
    void _vti_writer_micro_boolean(string outputname);
    void _vti_writer_atom (string outputname);
    void _vti_reader(string inputname);

    void _vtr_writer(string outputname);
    void _vtr_reader(string inputname);
    double meani(const int i);



cdef class pyState:
  cdef State* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in, m_in, evaluate_mean = None):
    # switch for evaluate_mean value
    if (evaluate_mean is None):
      self.thisptr = new State (deref(mesh_in.thisptr), deref(param_in.thisptr), ctypes.addressof(m_in.arr))  
    else:
      self.thisptr = new State (deref(mesh_in.thisptr), deref(param_in.thisptr), ctypes.addressof(m_in.arr), ctypes.addressof(evaluate_mean.arr))  
    #af.device.lock_array(m_in)#This does not avoid memory corruption caused by double free
  def __dealloc__(self): # causes segfault on every cleanup
    del self.thisptr
    self.thisptr = NULL
  @property
  def t(self):
    return self.thisptr.t
  @t.setter
  def t(self, value):
    self.thisptr.t = value
  def pythisptr(self):
      return <size_t><void*>self.thisptr

  def py_vti_writer_micro(self, outputname):
    self.thisptr._vti_writer_micro( outputname.encode('utf-8')) 
  def py_vti_writer_micro_boolean(self, outputname):
    self.thisptr._vti_writer_micro_boolean( outputname.encode('utf-8')) 
  def py_vti_writer_atom(self, outputname):
    self.thisptr._vti_writer_atom( outputname.encode('utf-8')) 
  def py_vti_reader(self, outputname):
    self.thisptr._vti_reader( outputname.encode('utf-8')) 

  def py_vtr_writer(self, outputname):
    self.thisptr._vtr_writer( outputname.encode('utf-8')) 
  def py_vtr_reader(self, outputname):
    self.thisptr._vtr_reader( outputname.encode('utf-8')) 
  def set_alpha(self,value):
    self.thisptr.param.alpha=value

  @property
  def m(self):
    m1=af.Array()
    m_addr = self.thisptr.get_m_addr()
    m1.arr=ctypes.c_void_p(m_addr)
    return m1
  @m.setter
  def m(self, m_in):
    self.thisptr.set_m(ctypes.addressof(m_in.arr))

  def meanxyz(self, i):
    return self.thisptr.meani(i)


cdef extern from "../../src/integrators/new_llg.hpp":
  cdef cppclass NewLlg:
    NewLlg (vector[shared_ptr[LLGTerm]] vector_in);
    vector[shared_ptr[LLGTerm]] llgterms;
    void step(State& state);
    double E(const State& state);
    long int get_fheff_addr(const State& state);
    #double cpu_time();
    #double h_stepped_;
    #State state0;

cdef class pyLLG:
  cdef NewLlg* thisptr
  def __cinit__(self, *args):
    cdef vector[shared_ptr[LLGTerm]] vector_in
    for arg in args:
      vector_in.push_back(shared_ptr[LLGTerm] (<LLGTerm*><size_t>arg.pythisptr()))
    self.thisptr = new NewLlg (vector_in)  
  # TODO leads to segfault on cleanup, compiler warning eleminated by adding virtual destructor in adaptive_rk.hpp
  # NOTE not happening in minimizer class as it is not derived (i guess)
  #def __dealloc__(self):
  #  del self.thisptr
  #  self.thisptr = NULL
  def llgstep(self, pyState state_in):
    self.thisptr.step(deref(state_in.thisptr))
  def get_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  #def print_stepsize(self):
  #  return self.thisptr.h_stepped_
  def get_fheff(self, pyState state):
    fheff_addr = self.thisptr.get_fheff_addr(deref(state.thisptr))
    fheff=af.Array()
    fheff.arr = ctypes.c_void_p(fheff_addr)
    return fheff
  #def cpu_time(self):
  #  return self.thisptr.cpu_time()
  #def set_state0_alpha(self,value):
  #  self.thisptr.state0.param.alpha=value
  def add_terms(self,*args):
    for arg in args:
      self.thisptr.llgterms.push_back(shared_ptr[LLGTerm] (<LLGTerm*><size_t>arg.pythisptr()))
    #cdef vector[shared_ptr[LLGTerm]] vector_in
    #for term in terms:
    #  vector_in.push_back(shared_ptr[LLGTerm] (<LLGTerm*><size_t>terms.pythisptr()))
    #self.thisptr = new NewLlg (vector_in)  
    
  

cdef extern from "../../src/llg_terms/micro_demag.hpp":
  cdef cppclass DemagSolver:
    DemagSolver (Mesh mesh_in, Param param_in);
    double E(const State& state);
    void print_Nfft();
    double get_cpu_time();

cdef class pyDemagSolver:
  cdef DemagSolver* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    self.thisptr = new DemagSolver (deref(mesh_in.thisptr), deref(param_in.thisptr))  
  #This would causes double free coruption!
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_Nfft(self):
    self.thisptr.print_Nfft()
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef extern from "../../src/llg_terms/micro_exch.hpp":
  cdef cppclass ExchSolver:
    ExchSolver (Mesh mesh_in, Param param_in);
    double E(const State& state);
    double get_cpu_time();

cdef class pyExchSolver:
  cdef ExchSolver* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    self.thisptr = new ExchSolver (deref(mesh_in.thisptr), deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef extern from "../../src/llg_terms/micro_anisotropy.hpp":
  cdef cppclass ANISOTROPY:
    ANISOTROPY(Mesh, Param);
    double E(const State& state);
    double get_cpu_time();

cdef class pyMicroAniso:
  cdef ANISOTROPY* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    self.thisptr = new ANISOTROPY (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr


cdef extern from "../../src/llg_terms/atomistic_demag.hpp":
  cdef cppclass ATOMISTIC_DEMAG:
    ATOMISTIC_DEMAG(Mesh);
    double E(const State& state);
    double get_cpu_time();

cdef class pyATOMISTIC_DEMAG:
  cdef ATOMISTIC_DEMAG* thisptr
  def __cinit__(self, pyMesh mesh_in):
    self.thisptr = new ATOMISTIC_DEMAG (deref(mesh_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef extern from "../../src/llg_terms/atomistic_anisotropy.hpp":
  cdef cppclass ATOMISTIC_ANISOTROPY:
    ATOMISTIC_ANISOTROPY(Mesh, Param);
    double E(const State& state);
    double get_cpu_time();

cdef class pyATOMISTIC_ANISOTROPY:
  cdef ATOMISTIC_ANISOTROPY* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    self.thisptr = new ATOMISTIC_ANISOTROPY (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef extern from "../../src/llg_terms/atomistic_exchange.hpp":
  cdef cppclass ATOMISTIC_EXCHANGE:
    ATOMISTIC_EXCHANGE(Mesh);
    double E(const State& state);
    double get_cpu_time();

cdef class AtomisticExchange:
  cdef ATOMISTIC_EXCHANGE* thisptr
  def __cinit__(self, pyMesh mesh_in):
    self.thisptr = new ATOMISTIC_EXCHANGE (deref(mesh_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef extern from "../../src/llg_terms/atomistic_dmi.hpp":
  cdef cppclass ATOMISTIC_DMI:
    ATOMISTIC_DMI(Mesh, Param);
    double E(const State& state);
    double get_cpu_time();

cdef class AtomisticDMI:
  cdef ATOMISTIC_DMI* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    self.thisptr = new ATOMISTIC_DMI (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef extern from "../../src/llg_terms/zee.hpp":
  cdef cppclass Zee:
    Zee (long int m_in);
    long int get_m_addr();
    double E(const State& state);
    double get_cpu_time();
    void set_xyz(const State&, const double x, const double y, const double z);

cdef class pyZee:
  cdef Zee* thisptr
  def __cinit__(self, array_in):
    self.thisptr = new Zee (ctypes.addressof(array_in.arr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,pyState state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def set_xyz(self,pyState state_in, x, y, z):
      self.thisptr.set_xyz(deref(state_in.thisptr), x, y, z)
  def pythisptr(self):
      return <size_t><void*>self.thisptr
  def get_zee(self):
    m1=af.Array()
    m_addr = self.thisptr.get_m_addr()
    m1.arr=ctypes.c_void_p(m_addr)
    return m1
#TODO 

cdef extern from "../../src/solvers/lbfgs_minimizer.hpp":
  cdef cppclass LBFGS_Minimizer:
    LBFGS_Minimizer(double tolerance_ , size_t maxIter_ , int verbose_ );
    vector[shared_ptr[LLGTerm]] llgterms_;
    #LBFGS_Minimizer(vector[shared_ptr[LLGTerm]] vector_in, double tolerance_ = 1e-6, size_t maxIter_ = 230, int verbose_ = 4);
    LBFGS_Minimizer(vector[shared_ptr[LLGTerm]] vector_in, double tolerance_, size_t maxIter_, int verbose_);
    double Minimize(State& state);
    double GetTimeCalcHeff();

cdef class pyLbfgsMinimizer:
  cdef LBFGS_Minimizer* thisptr
  def __cinit__(self, terms=[], tol = 1e-6, maxiter = 230):
    cdef vector[shared_ptr[LLGTerm]] vector_in
    if not terms:
      print("LBFGS_Minimizer: no terms provided, please add some either by providing a list terms=[...] or calling add_terms(*args)")
    else:
      for arg in terms:
        #print("Adding term", arg)
        vector_in.push_back(shared_ptr[LLGTerm] (<LLGTerm*><size_t>arg.pythisptr()))
      self.thisptr = new LBFGS_Minimizer (vector_in, tol, maxiter, 0) # TODO: WARNING: std::cout is not handled and leads to segfault!!! (setting verbose to 0 is a temporary fix) 
#TODO#  def __dealloc__(self): #causes segfault on GTO in cleanup
#TODO#    del self.thisptr
#TODO#    self.thisptr = NULL
  def add_terms(self,*args):
    for arg in args:
      self.thisptr.llgterms_.push_back(shared_ptr[LLGTerm] (<LLGTerm*><size_t>arg.pythisptr()))
  def delete_last_term(self):
    self.thisptr.llgterms_.pop_back()
  #not working as set_xyz is not pure virutal:# def set_zee_xyz(self, pyState state, i, x, y, z):
  #not working as set_xyz is not pure virutal:#     self.thisptr.llgterms_[i].set_xyz(deref(state.thisptr), x, y, z)
  #def __cinit__(self, tol = 1e-6, maxiter = 230, verbose = 4):
  #  self.thisptr = new LBFGS_Minimizer(tol, maxiter, verbose) # TODO handle default values 
  def pyMinimize(self, pyState state_in):
    return self.thisptr.Minimize(deref(state_in.thisptr))
  def pyGetTimeCalcHeff(self):
    return self.thisptr.GetTimeCalcHeff()


#cdef extern from "../../src/llg_terms/func.hpp":
#    double meani(const State& state_in, const int i);
#
#def pymeani(pyState state, int i):
#  return meani(deref(state.thisptr))

#cdef extern from "../../src/llg_terms/vtk_IO.hpp":
##  cdef void af_to_vti(array field, Mesh mesh, string outputname);
#  cdef void vti_writer_atom(array field, Mesh mesh, string outputname);
#def vti(field, mesh, outputname):
#  vti_writer_atom(field, mesh, outputname);

cdef extern from "../../src/param.hpp":
  cdef cppclass Param:
    Param();
    Param(double alpha, double T, double ms, double A, double D, double Ku1, double D_axis_x, double D_axis_y, double D_axis_z, double Ku1_axis_x, double Ku1_axis_y, double Ku1_axis_z, double p, double J_atom, double D_atom, double K_atom, double D_atom_axis_x , double D_atom_axis_y, double D_atom_axis_z, double K_atom_axis_x, double K_atom_axis_y, double K_atom_axis_z, bool hexagonal_close_packed, int mode, bool afsync);
    double mu0;
    double gamma;
    double alpha;
    double T;

    double ms;
    double A;
    double D;
    double Ku1;
    double D_axis[3];
    double Ku1_axis[3];

    double p;
    double J_atom;
    double D_atom;
    double K_atom;
    double D_atom_axis[3];
    double K_atom_axis[3];
    bool  hexagonal_close_packed;

    int mode;
    bool afsync;

cdef class pyParam:
  cdef Param* thisptr
  def __cinit__(self, alpha = 0., T = 0., ms = 0., A = 0., D = 0., Ku1 = 0., D_axis = [0.,0.,-1], Ku1_axis = [0.,0.,1.], p = 0., J_atom = 0., D_atom = 0., K_atom = 0., D_atom_axis = [0.,0.,1.], Ku1_atom_axis = [0.,0.,1.], bool hexagonal_close_packed = False, mode = 6, afsync = False):
    Ku1_axis_renormed = [x/(sqrt(Ku1_axis[0]**2 + Ku1_axis[1]**2 + Ku1_axis[2]**2)) for x in Ku1_axis]
    Ku1_atom_axis_renormed = [x/(sqrt(Ku1_atom_axis[0]**2 + Ku1_atom_axis[1]**2 + Ku1_atom_axis[2]**2)) for x in Ku1_atom_axis]
    D_axis_renormed = [x/(sqrt(D_axis[0]**2 + D_axis[1]**2 + D_axis[2]**2)) for x in D_axis]
    D_atom_axis_renormed = [x/(sqrt(D_atom_axis[0]**2 + D_atom_axis[1]**2 + D_atom_axis[2]**2)) for x in D_atom_axis]
    self.thisptr = new Param (alpha, T, ms, A, D, Ku1, D_axis_renormed[0], D_axis_renormed[1], D_axis_renormed[2], Ku1_axis_renormed[0], Ku1_axis_renormed[1], Ku1_axis_renormed[2], p, J_atom, D_atom, K_atom, D_atom_axis_renormed[0] , D_atom_axis_renormed[1], D_atom_axis_renormed[2], Ku1_atom_axis_renormed[0], Ku1_atom_axis_renormed[1], Ku1_atom_axis_renormed[2], hexagonal_close_packed , mode , afsync)
    #mu0 = self.thisptr.mu0
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL

  ## Getter and setter Functions
  @property # Note: constant by not providing setter
  def mu0(self):
    return self.thisptr.mu0

  @property
  def gamma(self):
    return self.thisptr.gamma

  @property
  def alpha(self):
    return self.thisptr.alpha
  @alpha.setter
  def alpha(self,value):
    self.thisptr.alpha=value

  @property
  def T(self):
    return self.thisptr.T
  @T.setter
  def T(self,value):
    self.thisptr.T=value

  # Micromagnetic
  @property
  def ms(self):
    return self.thisptr.ms
  @ms.setter
  def ms(self,value):
    self.thisptr.ms=value

  @property
  def A(self):
    return self.thisptr.A
  @A.setter
  def A(self,value):
    self.thisptr.A=value

  @property
  def D(self):
    return self.thisptr.D
  @D.setter
  def D(self,value):
    self.thisptr.D=value

  @property
  def D_axis(self):
    return self.thisptr.D_axis[0], self.thisptr.D_axis[1], self.thisptr.D_axis[2]
  @D_axis.setter
  def D_axis(self, values):
    self.thisptr.D_axis[0] = values[0]
    self.thisptr.D_axis[1] = values[1]
    self.thisptr.D_axis[2] = values[2]

  @property
  def Ku1(self):
    return self.thisptr.Ku1
  @Ku1.setter
  def Ku1(self,value):
    self.thisptr.Ku1=value

  @property
  def Ku1_axis(self):
    return self.thisptr.Ku1_axis[0], self.thisptr.Ku1_axis[1], self.thisptr.Ku1_axis[2]
  @Ku1_axis.setter
  def Ku1_axis(self, values):
    self.thisptr.Ku1_axis[0] = values[0]
    self.thisptr.Ku1_axis[1] = values[1]
    self.thisptr.Ku1_axis[2] = values[2]

  # Atomistic
  @property
  def p(self):
    return self.thisptr.p
  @p.setter
  def p(self,value):
    self.thisptr.p=value

  @property
  def J_atom(self):
    return self.thisptr.J_atom
  @J_atom.setter
  def J_atom(self,value):
    self.thisptr.J_atom=value

  @property
  def D_atom(self):
    return self.thisptr.D_atom
  @D_atom.setter
  def D_atom(self,value):
    self.thisptr.D_atom=value

  @property
  def Ku1_atom(self):
    return self.thisptr.K_atom
  @Ku1_atom.setter
  def Ku1_atom(self,value):
    self.thisptr.K_atom=value

  @property
  def D_atom_axis(self):
    return self.thisptr.D_atom_axis[0], self.thisptr.D_atom_axis[1], self.thisptr.D_atom_axis[2]
  @D_atom_axis.setter
  def D_atom_axis(self, values):
    self.thisptr.D_atom_axis[0] = values[0]
    self.thisptr.D_atom_axis[1] = values[1]
    self.thisptr.D_atom_axis[2] = values[2]

  @property
  def Ku1_atom_axis(self):
    return self.thisptr.K_atom_axis[0], self.thisptr.K_atom_axis[1], self.thisptr.K_atom_axis[2]
  @Ku1_atom_axis.setter
  def Ku1_atom_axis(self, values):
    self.thisptr.K_atom_axis[0] = values[0]
    self.thisptr.K_atom_axis[1] = values[1]
    self.thisptr.K_atom_axis[2] = values[2]
