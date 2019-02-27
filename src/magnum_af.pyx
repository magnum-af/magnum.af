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

from arrayfire import Array
from ctypes import addressof, c_void_p
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp cimport bool
from cython.operator cimport dereference as deref
from math import sqrt

from magnum_af_decl cimport ExchSolver as cExchSolver
from magnum_af_decl cimport Mesh as cMesh
from magnum_af_decl cimport Param as cParam
from magnum_af_decl cimport State as cState
from magnum_af_decl cimport ATOMISTIC_DMI
from magnum_af_decl cimport NewLlg
from magnum_af_decl cimport DemagSolver
from magnum_af_decl cimport ANISOTROPY
from magnum_af_decl cimport ATOMISTIC_DEMAG
from magnum_af_decl cimport ATOMISTIC_ANISOTROPY
from magnum_af_decl cimport ATOMISTIC_EXCHANGE
from magnum_af_decl cimport Zee
from magnum_af_decl cimport LBFGS_Minimizer
from magnum_af_decl cimport LLGTerm
from magnum_af_decl cimport array

#NOTE#@cython.embedsignature(True)# error: Cdef functions/classes cannot take arbitrary decorators. https://stackoverflow.com/questions/42668252/cython-cdef-class-not-displaying-doc-string-or-init-parameters
# Docstring does work, todo: check type etc. 
cdef class Mesh:
  """
  Class defining number of nodes and discretization.
  """
  def __init__(self, n0, n1, n2, dx, dy, dz): # todo: For dockstring, but currently does not change signature
    pass

  cdef cMesh* thisptr
  def __cinit__(self,int n0, int n1, int n2, double dx, double dy, double dz):
    self.thisptr = new cMesh(n0, n1, n2, dx, dy, dz)
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def n0(self):
    return self.thisptr.n0
  def print_n0(self):
    print (self.thisptr.n0)



cdef class State:
  cdef cState* thisptr
  def __cinit__(self, Mesh mesh_in, Param param_in, m_in, evaluate_mean = None):
    # switch for evaluate_mean value
    if (evaluate_mean is None):
      self.thisptr = new cState (deref(mesh_in.thisptr), deref(param_in.thisptr), addressof(m_in.arr))  
    else:
      self.thisptr = new cState (deref(mesh_in.thisptr), deref(param_in.thisptr), addressof(m_in.arr), addressof(evaluate_mean.arr))  
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
    m1=Array()
    m_addr = self.thisptr.get_m_addr()
    m1.arr=c_void_p(m_addr)
    return m1
  @m.setter
  def m(self, m_in):
    self.thisptr.set_m(addressof(m_in.arr))

  def meanxyz(self, i):
    return self.thisptr.meani(i)



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
  def llgstep(self, State state_in):
    self.thisptr.step(deref(state_in.thisptr))
  def get_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  #def print_stepsize(self):
  #  return self.thisptr.h_stepped_
  def get_fheff(self, State state):
    fheff_addr = self.thisptr.get_fheff_addr(deref(state.thisptr))
    fheff=Array()
    fheff.arr = c_void_p(fheff_addr)
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
    
  

cdef class pyDemagSolver:
  cdef DemagSolver* thisptr
  def __cinit__(self, Mesh mesh_in, Param param_in):
    self.thisptr = new DemagSolver (deref(mesh_in.thisptr), deref(param_in.thisptr))  
  #This would causes double free coruption!
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_Nfft(self):
    self.thisptr.print_Nfft()
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class ExchSolver:
  cdef cExchSolver* thisptr
  def __cinit__(self, Mesh mesh_in, Param param_in):
    self.thisptr = new cExchSolver (deref(mesh_in.thisptr), deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class pyMicroAniso:
  cdef ANISOTROPY* thisptr
  def __cinit__(self, Mesh mesh_in, Param param_in):
    self.thisptr = new ANISOTROPY (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr


cdef class pyATOMISTIC_DEMAG:
  cdef ATOMISTIC_DEMAG* thisptr
  def __cinit__(self, Mesh mesh_in):
    self.thisptr = new ATOMISTIC_DEMAG (deref(mesh_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class pyATOMISTIC_ANISOTROPY:
  cdef ATOMISTIC_ANISOTROPY* thisptr
  def __cinit__(self, Mesh mesh_in, Param param_in):
    self.thisptr = new ATOMISTIC_ANISOTROPY (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class AtomisticExchange:
  cdef ATOMISTIC_EXCHANGE* thisptr
  def __cinit__(self, Mesh mesh_in):
    self.thisptr = new ATOMISTIC_EXCHANGE (deref(mesh_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class AtomisticDMI:
  cdef ATOMISTIC_DMI* thisptr
  def __cinit__(self, Mesh mesh_in, Param param_in):
    self.thisptr = new ATOMISTIC_DMI (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class pyZee:
  cdef Zee* thisptr
  def __cinit__(self, array_in):
    self.thisptr = new Zee (addressof(array_in.arr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def set_xyz(self,State state_in, x, y, z):
      self.thisptr.set_xyz(deref(state_in.thisptr), x, y, z)
  def pythisptr(self):
      return <size_t><void*>self.thisptr
  def get_zee(self):
    m1=Array()
    m_addr = self.thisptr.get_m_addr()
    m1.arr=c_void_p(m_addr)
    return m1
#TODO 

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
  #not working as set_xyz is not pure virutal:# def set_zee_xyz(self, State state, i, x, y, z):
  #not working as set_xyz is not pure virutal:#     self.thisptr.llgterms_[i].set_xyz(deref(state.thisptr), x, y, z)
  #def __cinit__(self, tol = 1e-6, maxiter = 230, verbose = 4):
  #  self.thisptr = new LBFGS_Minimizer(tol, maxiter, verbose) # TODO handle default values 
  def pyMinimize(self, State state_in):
    return self.thisptr.Minimize(deref(state_in.thisptr))
  def pyGetTimeCalcHeff(self):
    return self.thisptr.GetTimeCalcHeff()

cdef class Param:
  cdef cParam* thisptr
  def __cinit__(self, alpha = 0., T = 0., ms = 0., A = 0., D = 0., Ku1 = 0., D_axis = [0.,0.,-1], Ku1_axis = [0.,0.,1.], p = 0., J_atom = 0., D_atom = 0., K_atom = 0., D_atom_axis = [0.,0.,1.], Ku1_atom_axis = [0.,0.,1.], bool hexagonal_close_packed = False, mode = 6, afsync = False):
    Ku1_axis_renormed = [x/(sqrt(Ku1_axis[0]**2 + Ku1_axis[1]**2 + Ku1_axis[2]**2)) for x in Ku1_axis]
    Ku1_atom_axis_renormed = [x/(sqrt(Ku1_atom_axis[0]**2 + Ku1_atom_axis[1]**2 + Ku1_atom_axis[2]**2)) for x in Ku1_atom_axis]
    D_axis_renormed = [x/(sqrt(D_axis[0]**2 + D_axis[1]**2 + D_axis[2]**2)) for x in D_axis]
    D_atom_axis_renormed = [x/(sqrt(D_atom_axis[0]**2 + D_atom_axis[1]**2 + D_atom_axis[2]**2)) for x in D_atom_axis]
    self.thisptr = new cParam (alpha, T, ms, A, D, Ku1, D_axis_renormed[0], D_axis_renormed[1], D_axis_renormed[2], Ku1_axis_renormed[0], Ku1_axis_renormed[1], Ku1_axis_renormed[2], p, J_atom, D_atom, K_atom, D_atom_axis_renormed[0] , D_atom_axis_renormed[1], D_atom_axis_renormed[2], Ku1_atom_axis_renormed[0], Ku1_atom_axis_renormed[1], Ku1_atom_axis_renormed[2], hexagonal_close_packed , mode , afsync)
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
