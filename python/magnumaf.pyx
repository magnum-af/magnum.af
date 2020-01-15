"""@package magnum.af
A finite differences GPU-accelerated micromagnetic and atomistic simulation software.
"""
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

import arrayfire as af

from ctypes import addressof, c_void_p
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp cimport bool
from cython.operator cimport dereference as deref
from math import sqrt
from math import pi
from numpy import zeros as np_zeros

from magnumaf_decl cimport Mesh as cMesh
from magnumaf_decl cimport NonequispacedMesh as cNonequispacedMesh
from magnumaf_decl cimport State as cState
from magnumaf_decl cimport Controller as cController
from magnumaf_decl cimport LLGIntegrator as cLLGIntegrator
from magnumaf_decl cimport DemagField as cDemagField
from magnumaf_decl cimport UniaxialAnisotropyField as cUniaxialAnisotropyField
from magnumaf_decl cimport NonequiUniaxialAnisotropyField as cNonequiUniaxialAnisotropyField
from magnumaf_decl cimport ExchangeField as cExchangeField
from magnumaf_decl cimport SparseExchangeField as cSparseExchangeField
from magnumaf_decl cimport NonequiExchangeField as cNonequiExchangeField
from magnumaf_decl cimport SpinTransferTorqueField as cSpinTransferTorqueField
from magnumaf_decl cimport RKKYExchangeField as cRKKYExchangeField
from magnumaf_decl cimport DmiField as cDmiField

from magnumaf_decl cimport AtomisticDipoleDipoleField as cAtomisticDipoleDipoleField
from magnumaf_decl cimport AtomisticExchangeField as cAtomisticExchangeField
from magnumaf_decl cimport AtomisticUniaxialAnisotropyField as cAtomisticUniaxialAnisotropyField
from magnumaf_decl cimport AtomisticDmiField as cAtomisticDmiField
from magnumaf_decl cimport ExternalField as cExternalField
from magnumaf_decl cimport AtomisticExternalField as cAtomisticExternalField
from magnumaf_decl cimport LBFGS_Minimizer as cLBFGS_Minimizer
from magnumaf_decl cimport LLGTerm as cLLGTerm

from magnumaf_decl cimport pywrap_vti_writer_micro as cpywrap_vti_writer_micro
from magnumaf_decl cimport String as cString

def array_from_addr(array_addr):
    array=af.Array()
    array.arr=c_void_p(array_addr)
    return array

class Util:
    @staticmethod
    def normalize(a):
        """Expects an vector array and returns the array normalized to 1."""
        norm_a = af.tile(af.sqrt(af.sum(a*a, 3)), 1, 1, 1, 3)
        normalized = a/norm_a
        af.replace(normalized, norm_a != 0, 0)
        return normalized

    @staticmethod
    def normed_homogeneous_field(nx = 1, ny = 1, nz = 1, axis=[1, 0, 0], factor = 1.):
        """Returns a homogeneous field of dimension [nx, ny, nz, 3] pointing into the direction of axis and normed to factor."""
        norm = sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
        array = af.constant(0.0, 1, 1, 1, 3, dtype=af.Dtype.f64)
        array [0, 0, 0, 0] = factor * axis[0]/norm
        array [0, 0, 0, 1] = factor * axis[1]/norm
        array [0, 0, 0, 2] = factor * axis[2]/norm
        return af.tile(array, nx, ny, nz)

    @staticmethod
    def disk(nx, ny, nz, axis=[1, 0, 0], return_ncells = False):
        norm = sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
        n_cells=0
        m = np_zeros((nx, ny, nz, 3))
        for ix in range (0, nx):
            for iy in range(0, ny):
                for iz in range(0, nz):
                    a= nx/2
                    b= ny/2
                    rx=ix-nx/2.
                    ry=iy-ny/2.
                    r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
                    if(r<=1):
                            m[ix, iy, iz, 0]=axis[0]/norm
                            m[ix, iy, iz, 1]=axis[1]/norm
                            m[ix, iy, iz, 2]=axis[2]/norm
                            n_cells = n_cells +1
        if return_ncells == True:
            return af.from_ndarray(m), n_cells
        else:
            return af.from_ndarray(m)

    @staticmethod
    def vortex(mesh, positive_z = True):
        """Returns a vortex configuration of dimension [mesh.nx, mesh.ny, mesh.nz, 3]. The option positive_z decides whether the vortex core points in positive (=True) or negative (=False) z-direction."""
        m = np_zeros((mesh.nx, mesh.ny, mesh.nz, 3));
        for ix in range (0, mesh.nx):
            for iy in range(0, mesh.ny):
                rx=float(ix)-mesh.nx/2.
                ry=float(iy)-mesh.ny/2.
                r = sqrt(pow(rx, 2)+pow(ry, 2))

                if r < mesh.nx/2.:
                    for iz in range(0, mesh.nz):
                        if r==0.:
                            if positive_z==True:
                                m[ix, iy, :, 2]= 1.
                            else:
                                m[ix, iy, :, 2]= -1
                        else:
                            m[ix, iy, :, 0]=-ry/r
                            m[ix, iy, :, 1]= rx/r
                            if positive_z==True:
                                m[ix, iy, :, 2]= sqrt(mesh.nx)/r
                            else:
                                m[ix, iy, :, 2]= - sqrt(mesh.nx)/r
                        norm = sqrt(m[ix, iy, iz, 0]**2+m[ix, iy, iz, 1]**2+m[ix, iy, iz, 2]**2)
                        m[ix, iy, iz, :]=m[ix, iy, iz, :]/norm
        return af.from_ndarray(m)

    @classmethod
    def sum_of_difference_of_abs(cls, a, b):
        return af.sum(af.sum(af.sum(af.sum(af.abs(a)-af.abs(b), 0), 1), 2), 3).scalar()

    @classmethod
    def test_sum_of_difference_of_abs(cls, a, b, verbose = True):
            c = cls.sum_of_difference_of_abs(a, b)
            if (c != 0.):
                    if (verbose == True):
                            print ("Error")
                    return False
            else:
                    if (verbose == True):
                            print ("Success")
                    return True

    @staticmethod
    def write_vti(afarray, dx, dy, dz, filename):
        cpywrap_vti_writer_micro(addressof(afarray.arr), dx, dy, dz, filename.encode('utf-8'))

# For adding methods as properties (e.g. the State class method m_partial as attribute)
# From http://code.activestate.com/recipes/440514-dictproperty-properties-for-dictionary-attributes/
class dictproperty(object):

    class _proxy(object):

        def __init__(self, obj, fget, fset, fdel):
            self._obj = obj
            self._fget = fget
            self._fset = fset
            self._fdel = fdel

        def __getitem__(self, key):
            if self._fget is None:
                raise TypeError, "can't read item"
            return self._fget(self._obj, key)

        def __setitem__(self, key, value):
            if self._fset is None:
                raise TypeError, "can't set item"
            self._fset(self._obj, key, value)

        def __delitem__(self, key):
            if self._fdel is None:
                raise TypeError, "can't delete item"
            self._fdel(self._obj, key)

    def __init__(self, fget=None, fset=None, fdel=None, doc=None):
        self._fget = fget
        self._fset = fset
        self._fdel = fdel
        self.__doc__ = doc

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return self._proxy(obj, self._fget, self._fset, self._fdel)


#NOTE#@cython.embedsignature(True)# error: Cdef functions/classes cannot take arbitrary decorators. https://stackoverflow.com/questions/42668252/cython-cdef-class-not-displaying-doc-string-or-init-parameters
# Docstring does work, todo: check type etc.
cdef class Mesh:
    """
    The Mesh object stores information about the discretization of our sample.


    Parameters
    ----------
    nx : int
        Number of cells in the x-direction
    ny : int
        Number of cells in the y-direction
    nz : int
        Number of cells in the z-direction
    dx : float
        Distance of one cell in the x-direction
    dy : float
        Distance of one cell in the y-direction
    dz : float
        Distance of one cell in the z-direction


    Attributes
    ----------
    All parameters are stored in equally named attributes


    Examples
    ----------
    mesh = Mesh(100, 25, 1, 1e-9, 1e-9, 1e-9)
    mesh = Mesh(nx = 100, ny = 25, nz = 1, dx = 1e-9, dy = 1e-9, dz = 1e-9)
    print(mesh.nx)
    """
    def __init__(self, nx, ny, nz, dx, dy, dz): # TODO: For dockstring, but currently does not change signature
        pass

    cdef cMesh* _thisptr
    cdef object owner # None if this is our own # From [1]
    def __cinit__(self, int nx, int ny, int nz, double dx, double dy, double dz):
        self._thisptr = new cMesh(nx, ny, nz, dx, dy, dz)
        owner = None # see [1]
    cdef set_ptr(self, cMesh* ptr, owner):
        if self.owner is None:
            del self._thisptr
        self._thisptr = ptr
        self.owner = owner
    def __dealloc__(self):
        if self.owner is None: # only free if we own it: see [1]
            del self._thisptr
            self._thisptr = NULL
    def print_nx(self):
        print (self._thisptr.n0)

    @property
    def nx(self):
        return self._thisptr.n0
    @nx.setter
    def nx(self, value):
        self._thisptr.n0=value

    @property
    def ny(self):
        return self._thisptr.n1
    @ny.setter
    def ny(self, value):
        self._thisptr.n1=value

    @property
    def nz(self):
        return self._thisptr.n2
    @nz.setter
    def nz(self, value):
        self._thisptr.n2=value

    @property
    def dx(self):
        return self._thisptr.dx
    @dx.setter
    def dx(self, value):
        self._thisptr.dx=value

    @property
    def dy(self):
        return self._thisptr.dy
    @dy.setter
    def dy(self, value):
        self._thisptr.dy=value

    @property
    def dz(self):
        return self._thisptr.dz
    @dz.setter
    def dz(self, value):
        self._thisptr.dz=value



cdef class NonequispacedMesh:
    """
    Nonequispaced Mesh object.
    """
    def __init__(self, nx, ny, dx, dy, z_spacing):
        pass

    cdef cNonequispacedMesh* _thisptr
    cdef object owner # None if this is our own # From [1]
    #def __cinit__(self, int nx, int ny, double dx, double dy, vector[double] z_spacing):
    def __cinit__(self, int nx, int ny, double dx, double dy, z_spacing):
        cdef vector[double] z_spacing_cvec
        for val in z_spacing:
            z_spacing_cvec.push_back(val)
        self._thisptr = new cNonequispacedMesh(nx, ny, dx, dy, z_spacing)
        owner = None # see [1]
    cdef set_ptr(self, cNonequispacedMesh* ptr, owner):
        if self.owner is None:
            del self._thisptr
        self._thisptr = ptr
        self.owner = owner
    def __dealloc__(self):
        if self.owner is None: # only free if we own it: see [1]
            del self._thisptr
            self._thisptr = NULL

    @property
    def nx(self):
        return self._thisptr.nx
    @nx.setter
    def nx(self, value):
        self._thisptr.nx=value

    @property
    def ny(self):
        return self._thisptr.ny
    @ny.setter
    def ny(self, value):
        self._thisptr.ny=value

    @property
    def nz(self):
        return self._thisptr.nz
    @nz.setter
    def nz(self, value):
        self._thisptr.nz=value

    @property
    def dx(self):
        return self._thisptr.dx
    @dx.setter
    def dx(self, value):
        self._thisptr.dx=value

    @property
    def dy(self):
        return self._thisptr.dy
    @dy.setter
    def dy(self, value):
        self._thisptr.dy=value

    @property
    def z_spacing(self):
        values = []
        for val in range(self._thisptr.nz):
            values.append(self._thisptr.z_spacing[val])
        return values
    #only add if needed, otherwise write only:
    #@z_spacing.setter
    #def z_spacing(self, values):
    #    cdef vector[double] cvec
    #    for val in values:
    #        cvec.push_back(val)
    #    self._thisptr.z_spacing = cvec


cdef class State:
    """
    The State object represents one realization of a magnetization configuration on a defined mesh.

    Parameters
    ----------
    mesh : Mesh
        The discretization is passed by a Mesh object
    Ms : float or af.array
            Saturation magnetization either set globally (float) or at each node (af.array)
    m : af.array
        Magnetization configuration which defines a vector at each node.
        The expected size is [mesh.nx, mesh.ny, mesh.nz, 3] and dtype=af.Dtype.f64
    verbose : bool(True)
        Print info messages, defaults to true
    mute_warning : bool(False)
        Mutes warnings, defaults to false

    Attributes
    ----------
    mesh : Mesh
    Ms : float
        Only set when input parameter Ms is float
    Ms_field : af.array
        Only set when input parameter Ms is af.array
    m : af.array [mesh.nx, mesh.ny, mesh.nz, 3, dtype=af.Dtype.f64]
    t : float
        Physical time in [s] manipulated by the LLGIntegrator object and initialized to 0

    Methods
    -------
    write_vti(outputname) : str
        Writes the current magnetization m into the file 'outputname.vti'
    normalize()
        Normalizes m to 1 for every node
    m_mean(i = None) : int
        Calculates the average magnetization along all three dimensions (i = None) or along dimension i = {0, 1, 2}

    Examples
    --------
    mesh = Mesh(1, 1, 1, 0.1, 0.1, 0.1)
    m0 = af.constant(0.0, nx, ny, nz, 3, dtype=af.Dtype.f64)
    m0[:, :, :, 0] = 1
    state = State(mesh, 8e5, m0)
    """
    cdef cState* _thisptr
    def __cinit__(self, Mesh mesh, Ms, m, verbose = True, mute_warning = False):
        # switch for evaluate_mean value
        if hasattr(Ms, 'arr'):
            self._thisptr = new cState (deref(mesh._thisptr), <long int> addressof(Ms.arr), <long int> addressof(m.arr), <bool> verbose, <bool> mute_warning)
        else:
            self._thisptr = new cState (deref(mesh._thisptr), <double> Ms, <long int> addressof(m.arr), <bool> verbose, <bool> mute_warning)
        #af.device.lock_array(m_in)#This does not avoid memory corruption caused by double free
    def __dealloc__(self): # causes segfault on every cleanup
        del self._thisptr
        self._thisptr = NULL
    @property
    def t(self):
        return self._thisptr.t
    @t.setter
    def t(self, value):
        self._thisptr.t = value
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    @property
    def Ms(self):
        return self._thisptr.Ms
    #@Ms.setter
    #def Ms(self, value):
    #  self._thisptr.Ms=value


    def write_vti(self, outputname):
        self._thisptr.write_vti( outputname.encode('utf-8'))
    def write_vti_boolean(self, outputname):
        self._thisptr._vti_writer_micro_boolean( outputname.encode('utf-8'))
    def write_vti_atomistic(self, outputname):
        self._thisptr._vti_writer_atom( outputname.encode('utf-8'))
    def read_vti(self, outputname):
        self._thisptr._vti_reader( outputname.encode('utf-8'))

    #def py_vtr_writer(self, outputname):
    #  self._thisptr._vtr_writer( outputname.encode('utf-8'))
    #def py_vtr_reader(self, outputname):
    #  self._thisptr._vtr_reader( outputname.encode('utf-8'))
    def normalize(self):
        self._thisptr.Normalize()

    property mesh:
        def __get__(self):
            mesh = Mesh(0, 0, 0, 0, 0, 0)
            mesh.set_ptr(&self._thisptr.mesh, self)
            return mesh
        def __set__(self, Mesh mesh):
            self._thisptr.mesh = deref(mesh._thisptr)

    property nonequimesh:
        def __set__(self, NonequispacedMesh ne_mesh):
            self._thisptr.nonequimesh = deref(ne_mesh._thisptr)

    @property
    def m(self):
        return array_from_addr(self._thisptr.get_m_addr())
    @m.setter
    def m(self, m_in):
        self._thisptr.set_m(addressof(m_in.arr))

    def get_m_partial(self, key):
        return array_from_addr(self._thisptr.get_m_addr())[key]

    def set_m_partial(self, key, value):
        if value.dims() == self.m[key].dims():
            temp = self.m # this copy is necessary, inline would lead to segfault
            temp[key] = value
            self._thisptr.set_m(addressof(temp.arr))
        else:
                print("Error: State.m_partial: Dimensions do not match. m_partial[key].dims()=", self.m[key].dims(), " != rhs.dims()=", value.dims(), ". Setting m_partial is ignored.")
    # setting dictionary as property
    m_partial = dictproperty(get_m_partial, set_m_partial, None)

    # Method for setting parts of m[key] to given af.array
    # Should be overloading m.setter, not completed
    #def __setitem__(self, key, m_in):
    #  temp = self.m
    #  temp[key] = m_in
    #  self._thisptr.set_m(addressof(temp.arr))
    #  print("min", self.m[key].dims(), "min", m_in.dims(), "m", self.m.dims(), "temp", temp.dims())
    #  #self._thisptr.set_m(addressof(m_in.arr))

    @property
    def Ms_field(self):
        return array_from_addr(self._thisptr.get_Ms_field())
    @Ms_field.setter
    def Ms_field(self, Ms_field):
        self._thisptr.set_Ms_field(addressof(Ms_field.arr))

    @property
    def steps(self):
        return self._thisptr.steps
    def m_mean(self, i = None):
        """
        Method calculating the average magnetization along all (i = None) or along a specific dimension ( i = {0, 1, 2})
        Parameter
        --------
        i : int (None)
        Returns
        ------
        <mx>, <my>, <mz>
            When i is omitted
        <mi>
            When i is either 0 = mx, 1 = my or 2 = mz
        """
        if(i == None):
            return self._thisptr.meani(0), self._thisptr.meani(1), self._thisptr.meani(2)
        else:
            return self._thisptr.meani(i)


cdef class LLGIntegrator:
    """
    LLGIntegrator(alpha, terms=[], mode="RKF45", hmin = 1e-15, hmax = 3.5e-10, atol = 1e-6, rtol = 1e-6)

    The LLGIntegrator object integrates a magnetization configuration according to the Landau–Lifshitz–Gilbert (LLG) equation

    Parameters
    ----------
    alpha : float
        The unitless damping constant in the LLG equation
    terms : [HeffTerm]
        A python list constisting of HeffTerm objects s.a. ExchangeField or DemagField
    mode : str
        Switch between Adapitve Runge Kutta schemes. Options are "RKF45", "DP45", "BS45", "DP78", "BS23"
    hmin : float
        Sets the minimum step size in seconds the integrator can take
    hmax : float
        Sets the minimum step size in seconds the integrator can take
    atol : float
        Sets the absolute tolerance for the adaptive Runge Kutta step-size
    rtol : float
        Sets the relative tolerance for the adaptive Runge Kutta step-size

    Attributes
    ----------
    alpha : float
        Stores the input parameter alpha

    Methods
    -------
    step(State)
        Make one Runge Kutta timestep integrating the LLG
        This manipulates state.m and state.t
    E(State) : float
        Calculates the micromagnetic energy of all terms for the magnetization state.m
    h(State) : af.array
        Returns the effective field H_eff for the magnetization state.m
    add_terms(*args)
        Adds an HeffTerm object (s.a. ExchangeField) to be included in the effective field
    relax(State, precision, ncalcE, nprint)
        Relaxes the magnetization until the energy difference between ncalcE steps is less than precision


    Examples
    ----------
    llg = LLGIntegrator (alpha = 0.02, terms = terms, mode = "RKF45")
    llg.step(state)
    print(llg.E(state))
    """
    cdef cLLGIntegrator* _thisptr
    def __cinit__(self, alpha, terms=[], mode="RKF45", hmin = 1e-15, hmax = 3.5e-10, atol = 1e-6, rtol = 1e-6):
        cdef vector[shared_ptr[cLLGTerm]] vector_in
        if not terms:
            print("LLGIntegrator: no terms provided, please add some either by providing a list LLGIntegrator(terms=[...]) or calling add_terms(*args) after declaration.")
        else:
            for arg in terms:
                vector_in.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg._get_thisptr()))
            self._thisptr = new cLLGIntegrator (alpha, vector_in, mode.encode('utf-8'), cController(hmin, hmax, atol, rtol))
    #def __dealloc__(self):
    #    # TODO maybe leads to segfault on cleanup, compiler warning eleminated by adding virtual destructor in adaptive_rk.hpp
    #    # NOTE is also problematic in minimizer class
    #    del self._thisptr
    #    self._thisptr = NULL
    def step(self, State state):
        self._thisptr.step(deref(state._thisptr))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def h(self, State state):
        return array_from_addr(self._thisptr.h_addr(deref(state._thisptr)))
    def add_terms(self, *args):
        for arg in args:
            self._thisptr.llgterms.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg._get_thisptr()))
    def relax(self, State state, precision = 1e-10, ncalcE = 100, nprint = 1000, verbose = True):
        """
        relax(State state, precision = 1e-10, ncalcE = 100, nprint = 1000)
            Relaxes the magnetization until the energy difference between ncalcE steps is less than precision
        """
        self._thisptr.relax(deref(state._thisptr), precision, ncalcE, nprint, verbose)
    @property
    def alpha(self):
        return self._thisptr.alpha
    @alpha.setter
    def alpha(self, value):
        self._thisptr.alpha=value

        #cdef vector[shared_ptr[cLLGTerm]] vector_in
        #for term in terms:
        #  vector_in.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>terms._get_thisptr()))
        #self._thisptr = new cLLGIntegrator (vector_in)

    #def print_stepsize(self):
    #  return self._thisptr.h_stepped_
    #def cpu_time(self):
    #  return self._thisptr.cpu_time()
    #def set_state0_alpha(self, value):
    #  self._thisptr.state0.material.alpha=value

cdef class String:
    """
    String method.
    """
    cdef cString* _thisptr
    #def __cinit__(self, State state, terms = [], n_interp=60, dt = 1e-13, LLGIntegrator llg):
    def __cinit__(self, State state, inputimages, n_interp, dt, LLGIntegrator llg):
        cdef vector[cState] vector_in
        if not inputimages:
            print("String: no States provided, please add at least two states for interpolation of initial path.")
        else:
            print(inputimages)
            for arg in inputimages:
                # swtich between state and array
                if hasattr(arg, 'arr'): # arg is expected to be array
                    mtemp = arg
                    vector_in.push_back(cState (cMesh(state.mesh.nx, state.mesh.ny, state.mesh.nz, state.mesh.dx, state.mesh.dy, state.mesh.dz), <double> state.Ms, <long int> addressof(mtemp.arr), <bool> True, <bool> True))
                else: # arg is expected to be State
                    mtemp = arg.m
                    if not arg.Ms_field.is_empty():
                        vector_in.push_back(cState (cMesh(state.mesh.nx, state.mesh.ny, state.mesh.nz, state.mesh.dx, state.mesh.dy, state.mesh.dz), <long int> addressof(arg.Ms_field.arr), <long int> addressof(mtemp.arr), <bool> True, <bool> True))
                    else:
                        vector_in.push_back(cState (cMesh(state.mesh.nx, state.mesh.ny, state.mesh.nz, state.mesh.dx, state.mesh.dy, state.mesh.dz), <double> state.Ms, <long int> addressof(mtemp.arr), <bool> True, <bool> True))

            self._thisptr = new cString (deref(state._thisptr), vector_in, n_interp, dt, deref(llg._thisptr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL

    def run(self, filepath, string_abort_rel_diff = 1e-12, string_abort_abs_diff = 1e-27, string_steps = 10000, verbose = True):
        return self._thisptr.run(filepath.encode('utf-8'), string_abort_rel_diff, string_abort_abs_diff, string_steps, <bool> verbose)

# base class for clearification (especially for heritage diagramms in docu)
cdef class HeffTerm:
    def h(self, State state):
        pass
    def E(self, State state):
        pass


cdef class DemagField(HeffTerm):
    """
    Demagnetization Field Object.

    mesh A mesh object
    verbose manage output
    caching enable reading and writing from cache
    nthreads set number of threads to use for the demag tensor computation. 0 uses all aviable threads
    """
    ## Demagnetization Field Object.
    #
    # @param mesh A mesh object
    # @param verbose manage output
    # @param caching enable reading and writing from cache
    # @param nthreads set number of threads to use for the demag tensor computation. 0 uses all aviable threads
    cdef cDemagField* _thisptr
    def __cinit__(self, Mesh mesh, verbose = False, caching = False, nthreads = 4):
        self._thisptr = new cDemagField (deref(mesh._thisptr), verbose, caching, nthreads)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def print_Nfft(self):
        self._thisptr.print_Nfft()
    ## Calculate energy contribution in [J]
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class ExchangeField(HeffTerm):
    cdef cExchangeField* _thisptr
    def __cinit__(self, A):
        if hasattr(A, 'arr'):
            self._thisptr = new cExchangeField (<long int> addressof(A.arr))
        else:
            self._thisptr = new cExchangeField (<double> A)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    ## Add when needed:
    # @property
    # def A(self):
    #       return self._thisptr.A
    # @A.setter
    # def A(self, value):
    #       self._thisptr.A=value
    # @property
    # def micro_A_field(self):
    #       return array_from_addr(self._thisptr.get_micro_A_field())
    # @micro_A_field.setter
    # def micro_A_field(self, micro_A_field_in):
    #       self._thisptr.set_micro_A_field(addressof(micro_A_field_in.arr))


cdef class SparseExchangeField(HeffTerm):
    cdef cSparseExchangeField* _thisptr
    def __cinit__(self, A, Mesh mesh, verbose = True):
        if hasattr(A, 'arr'):
            self._thisptr = new cSparseExchangeField (<long int> addressof(A.arr), deref(mesh._thisptr), <bool> verbose)
        else:
            self._thisptr = new cSparseExchangeField (<double> A, deref(mesh._thisptr), <bool> verbose)
            # Note: use <bool_t> instead of <bool> in case of ambiguous overloading error: https://stackoverflow.com/questions/29171087/cython-overloading-no-suitable-method-found
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class DmiField(HeffTerm):
    cdef cDmiField* _thisptr
    def __cinit__(self, D, D_axis = [0., 0., 1.]):
        if hasattr(D, 'arr'):
            self._thisptr = new cDmiField (<long int> addressof(D.arr), <double> D_axis[0], <double> D_axis[1], <double> D_axis[2])
        else:
            self._thisptr = new cDmiField (<double> D, <double> D_axis[0], <double> D_axis[1], <double> D_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class NonequiExchangeField(HeffTerm):
    cdef cNonequiExchangeField* _thisptr
    def __cinit__(self, A, NonequispacedMesh mesh, verbose = True):
        if hasattr(A, 'arr'):
            self._thisptr = new cNonequiExchangeField (<long int> addressof(A.arr), deref(mesh._thisptr), <bool> verbose)
        else:
            self._thisptr = new cNonequiExchangeField (<double> A, deref(mesh._thisptr), <bool> verbose)
            # Note: use <bool_t> instead of <bool> in case of ambiguous overloading error: https://stackoverflow.com/questions/29171087/cython-overloading-no-suitable-method-found
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr

cdef class UniaxialAnisotropyField(HeffTerm):
    cdef cUniaxialAnisotropyField* _thisptr
    def __cinit__(self, Ku1, Ku1_axis = [0, 0, 1]):
        if hasattr(Ku1, 'arr') and hasattr(Ku1_axis, 'arr'):
            self._thisptr = new cUniaxialAnisotropyField (<long int> addressof(Ku1.arr), <long int> addressof(Ku1_axis.arr))
        elif hasattr(Ku1, 'arr') :
            self._thisptr = new cUniaxialAnisotropyField (<long int> addressof(Ku1.arr), <double> Ku1_axis[0], <double> Ku1_axis[1], <double> Ku1_axis[2])
        elif hasattr(Ku1_axis, 'arr'):
            self._thisptr = new cUniaxialAnisotropyField (<double> Ku1, <long int> addressof(Ku1_axis.arr))
        else:
            self._thisptr = new cUniaxialAnisotropyField (<double> Ku1, <double> Ku1_axis[0], <double> Ku1_axis[1], <double> Ku1_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    @property
    def Ku1(self):
        return self._thisptr.Ku1
    #@Ku1.setter
    #def Ku1(self, value):
    #  self._thisptr.Ku1=value

    @property
    def Ku1_axis(self):
        return self._thisptr.get_ku1_axis(0), self._thisptr.get_ku1_axis(1), self._thisptr.get_ku1_axis(2)
    #@Ku1_axis.setter
    #def Ku1_axis(self, values):
    #  self._thisptr.Ku1_axis[0] = values[0]
    #  self._thisptr.Ku1_axis[1] = values[1]
    #  self._thisptr.Ku1_axis[2] = values[2]
    @property
    def Ku1_field(self):
        return array_from_addr(self._thisptr.get_Ku1_field())
    #@micro_Ku1_field.setter
    #def micro_Ku1_field(self, micro_Ku1_field_in):
    #  self._thisptr.set_micro_Ku1_field(addressof(micro_Ku1_field_in.arr))


cdef class NonequiUniaxialAnisotropyField(HeffTerm):
    cdef cNonequiUniaxialAnisotropyField* _thisptr
    def __cinit__(self, Ku1, Ku1_axis = [0, 0, 1]):
        if hasattr(Ku1, 'arr') and hasattr(Ku1_axis, 'arr'):
            self._thisptr = new cNonequiUniaxialAnisotropyField (<long int> addressof(Ku1.arr), <long int> addressof(Ku1_axis.arr))
        elif hasattr(Ku1, 'arr') :
            self._thisptr = new cNonequiUniaxialAnisotropyField (<long int> addressof(Ku1.arr), <double> Ku1_axis[0], <double> Ku1_axis[1], <double> Ku1_axis[2])
        elif hasattr(Ku1_axis, 'arr'):
            self._thisptr = new cNonequiUniaxialAnisotropyField (<double> Ku1, <long int> addressof(Ku1_axis.arr))
        else:
            self._thisptr = new cNonequiUniaxialAnisotropyField (<double> Ku1, <double> Ku1_axis[0], <double> Ku1_axis[1], <double> Ku1_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    @property
    def Ku1(self):
        return self._thisptr.Ku1
    #@Ku1.setter
    #def Ku1(self, value):
    #  self._thisptr.Ku1=value

    @property
    def Ku1_axis(self):
        return self._thisptr.get_ku1_axis(0), self._thisptr.get_ku1_axis(1), self._thisptr.get_ku1_axis(2)
    #@Ku1_axis.setter
    #def Ku1_axis(self, values):
    #  self._thisptr.Ku1_axis[0] = values[0]
    #  self._thisptr.Ku1_axis[1] = values[1]
    #  self._thisptr.Ku1_axis[2] = values[2]
    @property
    def Ku1_field(self):
        return array_from_addr(self._thisptr.get_Ku1_field())
    #@micro_Ku1_field.setter
    #def micro_Ku1_field(self, micro_Ku1_field_in):
    #  self._thisptr.set_micro_Ku1_field(addressof(micro_Ku1_field_in.arr))


cdef class AtomisticDipoleDipoleField(HeffTerm):
    cdef cAtomisticDipoleDipoleField* _thisptr
    def __cinit__(self, Mesh mesh):
        self._thisptr = new cAtomisticDipoleDipoleField (deref(mesh._thisptr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticUniaxialAnisotropyField(HeffTerm):
    cdef cAtomisticUniaxialAnisotropyField* _thisptr
    def __cinit__(self, double K_atom, K_atom_axis = [0., 0., 1.]):
    #def __cinit__(self, double K_atom, double K_atom_axis_x, double K_atom_axis_y, double K_atom_axis_z):
        self._thisptr = new cAtomisticUniaxialAnisotropyField (K_atom, K_atom_axis[0], K_atom_axis[1], K_atom_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticExchangeField(HeffTerm):
    cdef cAtomisticExchangeField* _thisptr
    def __cinit__(self, double J_atom):
        self._thisptr = new cAtomisticExchangeField (J_atom)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticDmiField(HeffTerm):
    cdef cAtomisticDmiField* _thisptr
    def __cinit__(self, double D_atom, D_atom_axis = [0., 0., -1.]):
        self._thisptr = new cAtomisticDmiField (D_atom, D_atom_axis[0], D_atom_axis[1], D_atom_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class ExternalField(HeffTerm):
    cdef cExternalField* _thisptr
    def __cinit__(self, array_in):
        self._thisptr = new cExternalField (addressof(array_in.arr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def set_homogeneous_field(self, x, y, z):
            self._thisptr.set_homogeneous_field(x, y, z)
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticExternalField(HeffTerm):
    cdef cAtomisticExternalField* _thisptr
    def __cinit__(self, array_in):
        self._thisptr = new cAtomisticExternalField (addressof(array_in.arr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.get_cpu_time()
    def set_homogeneous_field(self, x, y, z):
            self._thisptr.set_homogeneous_field(x, y, z)
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class SpinTransferTorqueField(HeffTerm):
    cdef cSpinTransferTorqueField* _thisptr
    def __cinit__(self, pol, nu_damp,  nu_field, j_e):
        self._thisptr = new cSpinTransferTorqueField (addressof(pol.arr), nu_damp, nu_field, j_e)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    @property
    def polarization_field(self):
        return array_from_addr(self._thisptr.polarization_field.get_array_addr())
    @polarization_field.setter
    def polarization_field(self, array):
        self._thisptr.polarization_field.set_array(addressof(array.arr))


cdef class RKKYExchangeField(HeffTerm):
    cdef cRKKYExchangeField* _thisptr
    def __cinit__(self, rkky_values, exchange_values, Mesh mesh, rkky_indices, verbose = True):
        self._thisptr = new cRKKYExchangeField (addressof(rkky_values.arr), addressof(exchange_values.arr), deref(mesh._thisptr), addressof(rkky_indices.arr), <bool> verbose)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def h(self, State state):
        return array_from_addr(self._thisptr.h_ptr(deref(state._thisptr)))
    def E(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class LBFGS_Minimizer:
    """
    The LBFGS_Minimizer object implements an Limited-memory Broyden–Fletcher–Goldfarb–Shanno (LBFGS) algorithm, minimizing a magentization configuratio with respect to the micromagnetic energy. This energy is obtained via the effective field given by the HeffTerm objects.

    Parameters
    ----------
    terms : [HeffTerm]
        List of HeffTerm objects
    tol : float (1e-6)
        Defines tolerance of the LBFGS algorithm
    maxiter : int (230)
        Defines maximum number of iterations in the LBFGS algorithm
    verbose : int(0)
        If > 0 enables output levels

    Methods
    -------
    minimize(State)
        Runns the minimization algorithm
    add_terms([HeffTerm])
        Adds a list of HeffTerm objects to the existing terms
    delete_last_term()
        Delets the lst HeffTerm object in the list

    Examples
    --------
    minimizer = Minimizer(terms, tol = 1e-6, maxiter = 230)
    minimier.minimize(state)
    """
    cdef cLBFGS_Minimizer* _thisptr
    def __cinit__(self, terms=[], tol = 1e-6, maxiter = 230, verbose = 0):
        cdef vector[shared_ptr[cLLGTerm]] vector_in
        if not terms:
            print("cLBFGS_Minimizer: no terms provided, please add some either by providing a list terms=[...] or calling add_terms(*args)")
        else:
            for arg in terms:
                vector_in.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg._get_thisptr()))
            self._thisptr = new cLBFGS_Minimizer (vector_in, tol, maxiter, verbose)
            # Note: std::cout is not handled and leads to segfault on some installs, hence c++ backend uses prinft

    #TODO
    #  def __dealloc__(self): #causes segfault on GTO in cleanup
    #      del self._thisptr
    #      self._thisptr = NULL
    def add_terms(self, *args):
        for arg in args:
            self._thisptr.llgterms_.push_back(shared_ptr[cLLGTerm] (<cLLGTerm*><size_t>arg._get_thisptr()))
    def delete_last_term(self):
        self._thisptr.llgterms_.pop_back()
    def minimize(self, State state):
        return self._thisptr.Minimize(deref(state._thisptr))
    def pyGetTimeCalcHeff(self):
        return self._thisptr.GetTimeCalcHeff()


class Constants:
    """
    Common physical constants. Values are obtained from CODATA/NIST.

    Attributes
    ----------
    mu0 : float
        [H/m] magnetic constant mu_0
    gamma : float
        [m A^-1 s^-1] gyromagnetic ratio gamma
    mu_b : float
        [J/T] Bohr magneton mu_bohr
    e : float
        [C] elementary charge e
    kb : float
        [J/K] Boltzmann constant kb
    hbar : float
        [J s] reduced Planck constant
    """
    mu0 = 4e-7 * pi

    gamma = 1.760859644e11 * mu0

    mu_b = 9.274009994e-24

    e = - 1.6021766208e-19

    kb = 1.38064852e-23

    hbar = 1.0545718e-34
