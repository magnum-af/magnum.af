#!python
#distutils: language = c++
#cython: language_level=3

# Links
# [1] https://stackoverflow.com/questions/33764094/cython-how-do-i-wrap-a-c-class-where-public-member-variables-are-custom-objec
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
from math import pi

from magnum_af_decl cimport Mesh as cMesh
from magnum_af_decl cimport Material as cParam
from magnum_af_decl cimport State as cState
from magnum_af_decl cimport LLGIntegrator as cLLGIntegrator
from magnum_af_decl cimport DemagField as cDemagField
from magnum_af_decl cimport UniaxialAnisotropyField as cUniaxialAnisotropyField
from magnum_af_decl cimport ExchangeField as cExchangeField
from magnum_af_decl cimport FieldlikeTorque as cFieldlikeTorque
from magnum_af_decl cimport DampinglikeTorque as cDampinglikeTorque
#TODO#from magnum_af_decl cimport DmiField as cDMI

from magnum_af_decl cimport AtomisticDipoleDipoleField as cAtomisticDipoleDipoleField
from magnum_af_decl cimport AtomisticExchangeField as cAtomisticExchangeField
from magnum_af_decl cimport AtomisticUniaxialAnisotropyField as cAtomisticUniaxialAnisotropyField
from magnum_af_decl cimport AtomisticDmiField as cAtomisticDmiField
from magnum_af_decl cimport Zee as cZee
from magnum_af_decl cimport LBFGS_Minimizer as cLBFGS_Minimizer
from magnum_af_decl cimport LLGTerm as cLLGTerm


#NOTE#@cython.embedsignature(True)# error: Cdef functions/classes cannot take arbitrary decorators. https://stackoverflow.com/questions/42668252/cython-cdef-class-not-displaying-doc-string-or-init-parameters
# Docstring does work, todo: check type etc. 
cdef class Mesh:
  """
  Class defining number of nodes and discretization.
  """
  def __init__(self, n0, n1, n2, dx, dy, dz): # todo: For dockstring, but currently does not change signature
    pass

  cdef cMesh* thisptr
  cdef object owner # None if this is our own # From [1]
  def __cinit__(self,int n0, int n1, int n2, double dx, double dy, double dz):
    self.thisptr = new cMesh(n0, n1, n2, dx, dy, dz)
    owner = None # see [1]
  cdef set_ptr(self, cMesh* ptr, owner):
    if self.owner is None:
      del self.thisptr
    self.thisptr = ptr
    self.owner = owner
  def __dealloc__(self):
    if self.owner is None: # only free if we own it: see [1]
      del self.thisptr
      self.thisptr = NULL
  def print_n0(self):
    print (self.thisptr.n0)

  @property
  def n0(self):
    return self.thisptr.n0
  @n0.setter
  def n0(self,value):
    self.thisptr.n0=value

  @property
  def n1(self):
    return self.thisptr.n1
  @n1.setter
  def n1(self,value):
    self.thisptr.n1=value

  @property
  def n2(self):
    return self.thisptr.n2
  @n2.setter
  def n2(self,value):
    self.thisptr.n2=value

  @property
  def dx(self):
    return self.thisptr.dx
  @dx.setter
  def dx(self,value):
    self.thisptr.dx=value

  @property
  def dy(self):
    return self.thisptr.dy
  @dy.setter
  def dy(self,value):
    self.thisptr.dy=value

  @property
  def dz(self):
    return self.thisptr.dz
  @dz.setter
  def dz(self,value):
    self.thisptr.dz=value



cdef class State:
  cdef cState* thisptr
  def __cinit__(self, Mesh mesh_in, Material param_in, m_in, evaluate_mean = None):
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
    self.thisptr.material.alpha=value
  def Normalize(self):
    self.thisptr.Normalize()

  property mesh:
    def __get__(self):
      mesh = Mesh(0,0,0,0,0,0)
      mesh.set_ptr(&self.thisptr.mesh,self)
      return mesh
    def __set__(self, Mesh mesh_in):
      self.thisptr.mesh = deref(mesh_in.thisptr)

  property material:
    def __get__(self):
      material = Material()
      material.set_ptr(&self.thisptr.material,self)
      return material
    def __set__(self, Material material_in):
      self.thisptr.material = deref(material_in.thisptr)

  @property
  def m(self):
    m1=Array()
    m_addr = self.thisptr.get_m_addr()
    m1.arr=c_void_p(m_addr)
    return m1
  @m.setter
  def m(self, m_in):
    self.thisptr.set_m(addressof(m_in.arr))

  @property
  def micro_Ms_field(self):
    micro_Ms_field1=Array()
    micro_Ms_field_addr = self.thisptr.get_micro_Ms_field()
    micro_Ms_field1.arr=c_void_p(micro_Ms_field_addr)
    return micro_Ms_field1
  @micro_Ms_field.setter
  def micro_Ms_field(self, micro_Ms_field_in):
    self.thisptr.set_micro_Ms_field(addressof(micro_Ms_field_in.arr))

  @property
  def micro_A_field(self):
    micro_A_field1=Array()
    micro_A_field_addr = self.thisptr.get_micro_A_field()
    micro_A_field1.arr=c_void_p(micro_A_field_addr)
    return micro_A_field1
  @micro_A_field.setter
  def micro_A_field(self, micro_A_field_in):
    self.thisptr.set_micro_A_field(addressof(micro_A_field_in.arr))

  @property
  def micro_Ku1_field(self):
    micro_Ku1_field1=Array()
    micro_Ku1_field_addr = self.thisptr.get_micro_Ku1_field()
    micro_Ku1_field1.arr=c_void_p(micro_Ku1_field_addr)
    return micro_Ku1_field1
  @micro_Ku1_field.setter
  def micro_Ku1_field(self, micro_Ku1_field_in):
    self.thisptr.set_micro_Ku1_field(addressof(micro_Ku1_field_in.arr))

  def meanxyz(self, i):
    return self.thisptr.meani(i)



cdef class LLGIntegrator:
  cdef cLLGIntegrator* thisptr
  def __cinit__(self, *args):
    cdef vector[shared_ptr[cLLGTerm]] vector_in
    for arg in args:
      vector_in.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg.pythisptr()))
    self.thisptr = new cLLGIntegrator (vector_in)  
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
  #  self.thisptr.state0.material.alpha=value
  def add_terms(self,*args):
    for arg in args:
      self.thisptr.llgterms.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg.pythisptr()))
    #cdef vector[shared_ptr[cLLGTerm]] vector_in
    #for term in terms:
    #  vector_in.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>terms.pythisptr()))
    #self.thisptr = new cLLGIntegrator (vector_in)  
    
  

cdef class DemagField:
  cdef cDemagField* thisptr
  def __cinit__(self, Mesh mesh_in, Material param_in):
    self.thisptr = new cDemagField (deref(mesh_in.thisptr), deref(param_in.thisptr))  
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

cdef class ExchangeField:
  cdef cExchangeField* thisptr
  def __cinit__(self, Mesh mesh_in, Material param_in):
    self.thisptr = new cExchangeField (deref(mesh_in.thisptr), deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class UniaxialAnisotropyField:
  cdef cUniaxialAnisotropyField* thisptr
  def __cinit__(self, Mesh mesh_in, Material param_in):
    self.thisptr = new cUniaxialAnisotropyField (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr


cdef class AtomisticDipoleDipoleField:
  cdef cAtomisticDipoleDipoleField* thisptr
  def __cinit__(self, Mesh mesh_in):
    self.thisptr = new cAtomisticDipoleDipoleField (deref(mesh_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class AtomisticUniaxialAnisotropyField:
  cdef cAtomisticUniaxialAnisotropyField* thisptr
  def __cinit__(self, Mesh mesh_in, Material param_in):
    self.thisptr = new cAtomisticUniaxialAnisotropyField (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class AtomisticExchangeField:
  cdef cAtomisticExchangeField* thisptr
  def __cinit__(self, Mesh mesh_in):
    self.thisptr = new cAtomisticExchangeField (deref(mesh_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class AtomisticDmiField:
  cdef cAtomisticDmiField* thisptr
  def __cinit__(self, Mesh mesh_in, Material param_in):
    self.thisptr = new cAtomisticDmiField (deref(mesh_in.thisptr),deref(param_in.thisptr))  
  def __dealloc__(self):
    del self.thisptr
    self.thisptr = NULL
  def print_E(self,State state_in):
    return self.thisptr.E(deref(state_in.thisptr))
  def cpu_time(self):
    return self.thisptr.get_cpu_time()
  def pythisptr(self):
      return <size_t><void*>self.thisptr

cdef class Zee:
  cdef cZee* thisptr
  def __cinit__(self, array_in):
    self.thisptr = new cZee (addressof(array_in.arr))  
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

cdef class LBFGS_Minimizer:
  cdef cLBFGS_Minimizer* thisptr
  def __cinit__(self, terms=[], tol = 1e-6, maxiter = 230):
    cdef vector[shared_ptr[cLLGTerm]] vector_in
    if not terms:
      print("cLBFGS_Minimizer: no terms provided, please add some either by providing a list terms=[...] or calling add_terms(*args)")
    else:
      for arg in terms:
        #print("Adding term", arg)
        vector_in.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg.pythisptr()))
      self.thisptr = new cLBFGS_Minimizer (vector_in, tol, maxiter, 0) # TODO: WARNING: std::cout is not handled and leads to segfault!!! (setting verbose to 0 is a temporary fix) 
#TODO#  def __dealloc__(self): #causes segfault on GTO in cleanup
#TODO#    del self.thisptr
#TODO#    self.thisptr = NULL
  def add_terms(self,*args):
    for arg in args:
      self.thisptr.llgterms_.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg.pythisptr()))
  def delete_last_term(self):
    self.thisptr.llgterms_.pop_back()
  #not working as set_xyz is not pure virutal:# def set_zee_xyz(self, State state, i, x, y, z):
  #not working as set_xyz is not pure virutal:#     self.thisptr.llgterms_[i].set_xyz(deref(state.thisptr), x, y, z)
  #def __cinit__(self, tol = 1e-6, maxiter = 230, verbose = 4):
  #  self.thisptr = new cLBFGS_Minimizer(tol, maxiter, verbose) # TODO handle default values 
  def pyMinimize(self, State state_in):
    return self.thisptr.Minimize(deref(state_in.thisptr))
  def pyGetTimeCalcHeff(self):
    return self.thisptr.GetTimeCalcHeff()

cdef class Material:
  cdef cParam* thisptr
  cdef object owner # None if this is our own # From [1]
  def __cinit__(self, alpha = 0., T = 0., ms = 0., A = 0., D = 0., Ku1 = 0., D_axis = [0.,0.,-1], Ku1_axis = [0.,0.,1.], p = 0., J_atom = 0., D_atom = 0., K_atom = 0., D_atom_axis = [0.,0.,1.], Ku1_atom_axis = [0.,0.,1.], bool hexagonal_close_packed = False, mode = 6, afsync = False):
    Ku1_axis_renormed = [x/(sqrt(Ku1_axis[0]**2 + Ku1_axis[1]**2 + Ku1_axis[2]**2)) for x in Ku1_axis]
    Ku1_atom_axis_renormed = [x/(sqrt(Ku1_atom_axis[0]**2 + Ku1_atom_axis[1]**2 + Ku1_atom_axis[2]**2)) for x in Ku1_atom_axis]
    D_axis_renormed = [x/(sqrt(D_axis[0]**2 + D_axis[1]**2 + D_axis[2]**2)) for x in D_axis]
    D_atom_axis_renormed = [x/(sqrt(D_atom_axis[0]**2 + D_atom_axis[1]**2 + D_atom_axis[2]**2)) for x in D_atom_axis]
    self.thisptr = new cParam (alpha, T, ms, A, D, Ku1, D_axis_renormed[0], D_axis_renormed[1], D_axis_renormed[2], Ku1_axis_renormed[0], Ku1_axis_renormed[1], Ku1_axis_renormed[2], p, J_atom, D_atom, K_atom, D_atom_axis_renormed[0] , D_atom_axis_renormed[1], D_atom_axis_renormed[2], Ku1_atom_axis_renormed[0], Ku1_atom_axis_renormed[1], Ku1_atom_axis_renormed[2], hexagonal_close_packed , mode , afsync)
    owner = None # see [1]
    #mu0 = self.thisptr.mu0
  cdef set_ptr(self, cParam* ptr, owner):
    if self.owner is None:
      del self.thisptr
    self.thisptr = ptr
    self.owner = owner
  def __dealloc__(self):
    if self.owner is None: # only free if we own it: see [1]
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

class Constants:
  """Common physical constants. Values were obtained from CODATA/NIST.
     Attributes:
     mu0 [H/m] magnetic constant mu_0
     gamma [m A^-1 s^-1] gyromagnetic ratio gamma
     mu_b [J/T] Bohr magneton mu_bohr
     e [C] elementary charge e
     kb [J/K] Boltzmann constant kb
     hbar [J s] reduced Planck constant"""
  mu0 = 4e-7 * pi

  gamma = 1.760859644e11 * mu0

  mu_b = 9.274009994e-24

  e = - 1.6021766208e-19

  kb = 1.38064852e-23

  hbar = 1.0545718e-34

cdef class FieldlikeTorque:
  cdef cFieldlikeTorque* thisptr
  def __cinit__(self, polarization_field, nu_field, j_e):
    self.thisptr = new cFieldlikeTorque (addressof(polarization_field.arr), nu_field, j_e) 

  @property
  def polarization_field(self):
    polarization_field1=Array()
    polarization_field_addr = self.thisptr.get_polarization_field_addr()
    polarization_field1.arr=c_void_p(polarization_field_addr)
    return polarization_field1
  @polarization_field.setter
  def polarization_field(self, polarization_field_in):
    self.thisptr.set_polarization_field(addressof(polarization_field_in.arr))
