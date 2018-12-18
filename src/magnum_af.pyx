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
from cython.operator cimport dereference as deref
import ctypes


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

cdef class pyMesh:
  cdef Mesh* thisptr
  def __cinit__(self,int a, int b, int c, double d, double e, double f):
    self.thisptr = new Mesh(a,b,c,d,e,f)
  def __dealloc__(self):
    del self.thisptr
  def n0(self):
    return self.thisptr.n0
  def print_n0(self):
    print (self.thisptr.n0)


cdef extern from "../../src/state.hpp":
  cdef cppclass State:
    State (Mesh mesh_in, Param param_in, long int m_in);
    State (Mesh mesh_in, Param param_in, long int aptr, long int evaluate_mean_ptr);
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
  def t(self):
    return self.thisptr.t
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

  def get_m(self):
    #af.set_backend("cpu")
    #af.set_device(0)
    #m1=af.Array()
    #m.arr = copy.deepcopy ()
    #m.arr = ctypes.c_void_p(m_addr)
    #af.sync()
    #af.device.lock_array(m)
    #af.device.lock_array(m2)
    m1=af.Array()
    m_addr = self.thisptr.get_m_addr()
    m1.arr=ctypes.c_void_p(m_addr)
    #print "test", m1.arr
    return m1
    #return m.copy()
    ##Alternative
    #return self.thisptr.get_m_addr()
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
  # TODO dealloc?
  #def __dealloc__(self):
  #  del self.thisptr
  def llgstep(self, pyState state_in):
    self.thisptr.step(deref(state_in.thisptr))
  def print_E(self,pyState state_in):
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

cdef extern from "../../src/llg_terms/micro_exch.hpp":
  cdef cppclass ExchSolver:
    ExchSolver (Mesh mesh_in, Param param_in);
    double E(const State& state);
    double get_cpu_time();

cdef class pyExchSolver:
  cdef ExchSolver* thisptr
  def __cinit__(self, pyMesh mesh_in, pyParam param_in):
    self.thisptr = new ExchSolver (deref(mesh_in.thisptr), deref(param_in.thisptr))  
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
  def __dealloc__(self):
    del self.thisptr
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
    double mu0;
    double gamma;
    double ms,A,alpha;
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
  #def __dealloc__(self):
  #  del self.thisptr

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
  def D_axis_x(self,value):
    self.thisptr.D_axis[0]=value
  def D_axis_y(self,value):
    self.thisptr.D_axis[1]=value
  def D_axis_z(self,value):
    self.thisptr.D_axis[2]=value
  def D_axis(self, *args):
    i=0
    for arg in args:
      self.thisptr.D_axis[i]=arg
      i=i+1
  def Ku1(self,value):
    self.thisptr.Ku1=value
  def Ku1_axis_x(self,value):
    self.thisptr.Ku1_axis[0]=value
  def Ku1_axis_y(self,value):
    self.thisptr.Ku1_axis[1]=value
  def Ku1_axis_z(self,value):
    self.thisptr.Ku1_axis[2]=value
  def Ku1_axis(self, *args):
    i=0
    for arg in args:
      self.thisptr.Ku1_axis[i]=arg
      i=i+1
  def J_atom(self,value):
    self.thisptr.J_atom=value
  def D_atom(self,value):
    self.thisptr.D_atom=value
  def D_atom_axis_x(self,value):
    self.thisptr.D_atom_axis[0]=value
  def D_atom_axis_y(self,value):
    self.thisptr.D_atom_axis[1]=value
  def D_atom_axis_z(self,value):
    self.thisptr.D_atom_axis[2]=value
  def D_atom_axis(self, *args):
    i=0
    for arg in args:
      self.thisptr.D_atom_axis[i]=arg
      i=i+1
  def K_atom(self,value):
    self.thisptr.K_atom=value
  def K_atom_axis_x(self,value):
    self.thisptr.K_atom_axis[0]=value
  def K_atom_axis_y(self,value):
    self.thisptr.K_atom_axis[1]=value
  def K_atom_axis_z(self,value):
    self.thisptr.K_atom_axis[2]=value
  def K_atom_axis(self, *args):
    i=0
    for arg in args:
      self.thisptr.K_atom_axis[i]=arg
      i=i+1

  def print_Ku1_axis(self):
    return self.thisptr.Ku1_axis[0], self.thisptr.Ku1_axis[1], self.thisptr.Ku1_axis[2]
  def print_Ku1_axis_x(self):
    return self.thisptr.Ku1_axis[0]
  def print_Ku1_axis_y(self):
    return self.thisptr.Ku1_axis[1]
  def print_Ku1_axis_z(self):
    return self.thisptr.Ku1_axis[2]

  def print_K_atom_axis_x(self):
    return self.thisptr.K_atom_axis[0]
  def print_K_atom_axis_y(self):
    return self.thisptr.K_atom_axis[1]
  def print_K_atom_axis_z(self):
    return self.thisptr.K_atom_axis[2]
  #Read only:
  def print_mu0(self):
    return self.thisptr.mu0
  def print_gamma(self):
    return self.thisptr.gamma
  def print_ms(self):
    return self.thisptr.ms
  def print_alpha(self):
    return self.thisptr.alpha
  def print_A(self):
    return self.thisptr.A
  def print_p(self):
    return self.thisptr.p
  def print_mode(self):
    return self.thisptr.mode
  def print_D(self):
    return self.thisptr.D
  def print_Ku1(self):
    return self.thisptr.Ku1
  def print_J_atom(self):
    return self.thisptr.J_atom
  def print_D_atom(self):
    return self.thisptr.D_atom
  def print_K_atom(self):
    return self.thisptr.K_atom
