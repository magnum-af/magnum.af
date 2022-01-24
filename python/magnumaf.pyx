# For Doxygen (package name must match filename):
## @package magnumafpy
# magnum.af python bindings.

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
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector
from libcpp cimport bool
from cython.operator cimport dereference as deref
from math import sqrt
from math import pi
import numpy as np
# from numpy import zeros as np_zeros

#from magnumafpy_decl cimport array_d3 as carray_d3
from magnumafpy_decl cimport Mesh as cMesh
from magnumafpy_decl cimport State as cState
from magnumafpy_decl cimport Controller as cController
from magnumafpy_decl cimport LLGIntegrator as cLLGIntegrator
from magnumafpy_decl cimport Stochastic_LLG as cStochastic_LLG
from magnumafpy_decl cimport LBFGS_Minimizer as cLBFGS_Minimizer
from magnumafpy_decl cimport StringMethod as cString

## FieldTerms
from magnumafpy_decl cimport FieldTerm as cFieldterm
# Micro
from magnumafpy_decl cimport ExternalField as cExternalField
from magnumafpy_decl cimport DemagField as cDemagField
from magnumafpy_decl cimport DemagFieldPBC as cDemagFieldPBC
from magnumafpy_decl cimport UniaxialAnisotropyField as cUniaxialAnisotropyField
from magnumafpy_decl cimport CubicAnisotropyField as cCubicAnisotropyField
from magnumafpy_decl cimport ExchangeField as cExchangeField
from magnumafpy_decl cimport SparseExchangeField as cSparseExchangeField
from magnumafpy_decl cimport ExchangeFieldPBC as cExchangeFieldPBC
from magnumafpy_decl cimport SpinTransferTorqueField as cSpinTransferTorqueField
from magnumafpy_decl cimport RKKYExchangeField as cRKKYExchangeField
from magnumafpy_decl cimport DmiField as cDmiField
from magnumafpy_decl cimport BulkDMIExchangeField as cBulkDMIExchangeField
from magnumafpy_decl cimport DMI_D2d_Field as cDMI_D2d_Field

# Nonequi
from magnumafpy_decl cimport NonequiMesh as cNonequiMesh
from magnumafpy_decl cimport NonequiDemagField as cNonequiDemagField
from magnumafpy_decl cimport NonequiExternalField as cNonequiExternalField
from magnumafpy_decl cimport NonequiUniaxialAnisotropyField as cNonequiUniaxialAnisotropyField
from magnumafpy_decl cimport NonequiExchangeField as cNonequiExchangeField

# Atomistic
from magnumafpy_decl cimport AtomisticDipoleDipoleField as cAtomisticDipoleDipoleField
from magnumafpy_decl cimport AtomisticExchangeField as cAtomisticExchangeField
from magnumafpy_decl cimport AtomisticUniaxialAnisotropyField as cAtomisticUniaxialAnisotropyField
from magnumafpy_decl cimport AtomisticDmiField as cAtomisticDmiField
from magnumafpy_decl cimport AtomisticExternalField as cAtomisticExternalField

# Util
from magnumafpy_decl cimport pywrap_vti_writer_micro as cpywrap_vti_writer_micro
from magnumafpy_decl cimport pywrap_vtr_writer as cpywrap_vtr_writer
from magnumafpy_decl cimport double_array3
from magnumafpy_decl cimport spacial_mean_in_region as cspacial_mean_in_region

# wrapping std::move
# runs with cython 0.29.14 (ubuntu 20.04 distro version), not working with cython 0.26 (ubuntu 18.04)
# with cython 0.29.21+ we could even only use:
# from libcpp.utility cimport move
# we move constructed std::vector<unique_ptr<Base>> instances into the wrapped ctors
# not moving would cause "error: static assertion failed: result type must be constructible from value type of input range"
# from https://github.com/cython/cython/issues/2169
cdef extern from * namespace "polyfill":
    """
    namespace polyfill {

    template <typename T>
    inline typename std::remove_reference<T>::type&& move(T& t) {
        return std::move(t);
    }

    template <typename T>
    inline typename std::remove_reference<T>::type&& move(T&& t) {
        return std::move(t);
    }

    }  // namespace polyfill
    """
    cdef T move[T](T)

# argparse
import argparse
import os
import sys
from shutil import copy

def parse():
    """
    Invokes argument parser.
    """
    parser = argparse.ArgumentParser(description='magnum.af simulation script.')
    parser.add_argument(
            '-o',
            '--outdir',
            type=str,
            default='output_' + os.path.basename(os.path.splitext(sys.argv[0])[0]) + '/',
            help="Output directory. Will be created and is accessible via 'parse.outdir'. Defaults to 'output_<scriptname>'.",
            )
    parser.add_argument(
            '-b',
            '--backend',
            type=str,
            help="'cuda', 'opencl' or 'cpu'. Select arrayfire backend via 'af.set_backend(b)'."
            )
    parser.add_argument(
            '-d',
            '--device',
            type=int,
            help="Set 'af.set_device(d)', e.g. used for selecting a GPU."
            )
    parser.add_argument(
            '-i',
            '--info',
            help='Prints af.info() to display active backend and device.',
            action='store_true' # defaults to 'False'
            )
    parser.add_argument(
            '-n',
            '--no-overwrite',
            help='Abort if outdir already exists. Prevents file overwriting.',
            action='store_true' # defaults to 'False'
            )
    parser.add_argument(
            '-N',
            '--nocopy',
            help='Skip copying this script into outdir/ as backup.',
            action='store_true' # defaults to 'False'
            )
    parser.add_argument(
            '-v',
            '--verbose',
            help="Make this parser verbose, printing parsing steps.",
            action='store_true'
            )
    parser.add_argument(
            'posargs',
            nargs=argparse.REMAINDER,
            help="Positional arguments, access via 'parser.posargs'.",
            )
    args = parser.parse_args()

    def print_if(switch, *args, **kwargs):
        """Prints preceding arguments if first argument evaluates to true"""
        if switch:
            print(*args, **kwargs)

    if not args.outdir.endswith('/'):
        args.outdir = args.outdir + '/'

    print_if(args.verbose, '--verbose: parse()=', args)

    if os.path.exists(args.outdir):
        if args.no_overwrite:
            raise RuntimeError("parse(): --no-overwrite: Output directory exists, aborting! Disable -n to overwrite.")
        else:
            print_if(args.verbose, "--outdir: Writing into existing directory '", args.outdir, "'")
    else:
        print_if(args.verbose, "--outdir: Dir does not exist, creating directory '", args.outdir, "'")
        os.makedirs(args.outdir, exist_ok = False)

    if args.backend is not None:
        print_if(args.verbose, "--backend: setting backend to", args.backend)
        af.set_backend(args.backend)

    if args.device is not None:
        print_if(args.verbose, "--device: setting device to", args.device)
        af.set_device(args.device)

    if args.nocopy is False:
        targetname = args.outdir + os.path.basename(sys.argv[0]) + '.copy'
        print_if(args.verbose, '--nocopy (not set): copying script "', os.path.basename(sys.argv[0]), '" into "', targetname, '"')
        copy(sys.argv[0], targetname)

    if args.info is True:
        af.info()

    return args

def parse_filepath():
    """
    Parses input and returns the parsed filepath using parser.pars_args().outdir
    """
    return parse().outdir

def array_from_addr(array_addr):
    array=af.Array()
    array.arr=c_void_p(array_addr)
    return array

def get_dim4(array):
    """Retruns list containing array.dims() as [n0, n1, n2, n3] """
    dim4 = [1] * 4
    for i,val in enumerate(array.dims(), start = 0):
        dim4[i] = val
    return dim4

def lookup_4d(values : af.Array, index : af.Array):
    """
    Performs af.lookup(values, af.flat(index)) and tiles to index.dims()
    values : [n_regions, 1, 1, 1]
    index  : [n0, n1, n2, n3]
    """
    look = af.lookup(values, af.flat(index))
    dims = get_dim4(index)
    return af.moddims(look, d0 = dims[0], d1 = dims[1], d2 = dims[2], d3 = dims[3])

def lookup(values : [float], index : af.Array, as_type = af.Dtype.f64):
    """
    Create array from lookup-table, unsing index as keys.
        values: lookup-table. Keys encoded as position.
                values[0] has key 0, values[1] has key 1, ...
        index:  Array holding keys to values.
                0 maps to values[0], 1 to values[1], ...
        returns: array containing values for respective keys.
    """
    values_array = af.array.Array(values).as_type(as_type)
    return lookup_4d(values_array, index)

def make_nonzeros_ones(array):
    """
    Zero valus stay zero.
    Non-zero values become Ones.
    Aternative names: isnotzero, make_binary, one_if_nonzero_else_zero
    """
    return af.iszero(af.iszero(array))


class Conversion:
    @staticmethod
    def Apm_to_Tesla(value_Apm):
        return value_Apm * Constants.mu0
    @staticmethod
    def Tesla_to_Apm(value_Tesla):
        return value_Tesla / Constants.mu0
    @staticmethod
    def J_to_eV(value_Joule):
        return value_Joule / Constants.e_abs
    @staticmethod
    def eV_to_J(value_eV):
        return value_eV * Constants.e_abs

class Util:
    @staticmethod
    def np_to_af(a):
        """
        Convert a numpy array to an arrayfire array.
        """
        return af.from_ndarray(a)
    @staticmethod
    def af_to_np(a):
        """
        Convert an arrayfire array to a numpy array.
        """
        return a.to_ndarray()
    @staticmethod
    def plot(plotfile = 'plot.gpi',
            outputfile = 'm.pdf',
            datafile = 'm.dat',
            xlabel = "t [s]",
            ylabel = '<m>',
            lines = ['u 1:2 w l t "m"'],
            gnuplot = True,
            evince = True,
            outputdir = None
            ):
        if outputdir is not None:
            plotfile = outputdir + plotfile
            outputfile = outputdir + outputfile
            datafile = outputdir + datafile
        f = open(plotfile, 'w')
        f.write('set terminal pdf\n')
        f.write('set output "' + outputfile + '"\n')
        f.write('set xlabel "' + xlabel + '"\n')
        f.write('set ylabel "' + ylabel + '"\n')
        f.write('plot ')
        for i, line in enumerate(lines):
            if i > 0:
                f.write(', ')
            f.write('"' + datafile + '" ' + line)
        f.write('\n')
        f.close()

        from os import system
        if gnuplot:
            system('gnuplot ' + plotfile)
        if evince:
            system('evince ' + outputfile)

    @staticmethod
    def spacial_mean(field):
        """Spacial mean of either a scalar or vector field"""
        mean = af.mean(af.mean(af.mean(field, dim=0), dim=1), dim=2)
        if len(mean.shape) <= 3:
            return mean.scalar()
        elif mean.shape[3] == 3:
            return mean[0, 0, 0, 0].scalar(), mean[0, 0, 0, 1].scalar(), mean[0, 0, 0, 2].scalar()
        else:
            raise RuntimeError('spacial_mean: input array has wrong dimension')

    @staticmethod
    def spacial_mean_in_region(vectorfield, region):
        """
        Evaluates mean values of vectorfield in every cell where region is not null.
        vectorfield: afarray [nx ny nz 3]
        region: afarray [nx ny nz 1]
            Each cell where region is not null is counted.
            Cells with zero are ignored.
        """
        a = cspacial_mean_in_region(addressof(vectorfield.arr), addressof(region.arr))
        return a[0], a[1], a[2]

    @staticmethod
    def normalize(a):
        """Expects an vector array and returns the array normalized to 1."""
        norm_a = af.tile(af.sqrt(af.sum(a*a, 3)), 1, 1, 1, 3)
        normalized = a/norm_a
        af.replace(normalized, norm_a != 0, 0)
        return normalized

    @staticmethod
    def disk(nx, ny, nz, axis=[1, 0, 0], return_ncells = False):
        norm = sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
        n_cells=0
        m = np.zeros((nx, ny, nz, 3))
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
        m = np.zeros((mesh.nx, mesh.ny, mesh.nz, 3));
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

    @staticmethod
    def sum_of_difference_of_abs(a, b):
        return af.sum(af.sum(af.sum(af.sum(af.abs(a)-af.abs(b), 0), 1), 2), 3).scalar()

    @staticmethod
    def test_sum_of_difference_of_abs(a, b, verbose = True):
            c = Util.sum_of_difference_of_abs(a, b)
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

    @staticmethod
    def write_vtr(afarray, NonequiMesh nonequi_mesh, filename):
        cpywrap_vtr_writer(addressof(afarray.arr), deref(nonequi_mesh._thisptr), filename.encode('utf-8'))

    @staticmethod
    def gto_gpu_renumeration(gpu_number):
        """Re-interprets gpu enumeration to be consistent with output of 'nvidia-smi' command. Returns abs(gpu_number -3) when active backend is 'cuda', else returns 0 for compatibility with cpu and opencl."""
        if af.library.get_active_backend() == 'cuda':
            return abs(gpu_number - 3)
        else:
            return 0

class Geometry:
    @staticmethod
    def xy_ellipse(nx, ny, nz, make_3d=False):
        """Geometry array representing an ellipse in xy plane with 1, else 0."""
        mesh = np.zeros((nx, ny, nz), dtype = np.bool)
        for ix in range (0, nx):
            for iy in range(0, ny):
                for iz in range(0, nz):
                    a= nx/2.
                    b= ny/2.
                    rx=ix-nx/2.
                    ry=iy-ny/2.
                    r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2)
                    if(r<=1):
                        mesh[ix, iy, iz]=1
        if make_3d:
            return af.tile(af.from_ndarray(mesh), 1, 1, 1, 3)
        else:
            return af.from_ndarray(mesh)

class Magnetization:
    """Functions returning normalized magnetizations. Can easily be combined with Geometry functions."""

    @staticmethod
    def isotropic(nx: int, ny: int, nz: int, dtype = af.Dtype.f64):
        """Returns an isotropic random distribution of unit vectors."""
        # random normal distribution of coordinates gives isotropic distribution of directions:
        return Util.normalize(af.randn(nx, ny, nz, 3, dtype))

    @staticmethod
    def homogeneous(nx: int, ny: int, nz: int, axis = [1, 0, 0], dtype=af.Dtype.f64):
        """Returns a normalized homogeneous field of dimension [nx, ny, nz, 3] pointing into the direction of axis."""
        norm = sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
        array = af.constant(0.0, 1, 1, 1, 3, dtype)
        array [0, 0, 0, 0] = axis[0]/norm
        array [0, 0, 0, 1] = axis[1]/norm
        array [0, 0, 0, 2] = axis[2]/norm
        return af.tile(array, nx, ny, nz)

class Field:
    @staticmethod
    def homogeneous(Hext, nx: int, ny: int, nz: int, dtype=af.Dtype.f64):
        """Returns a homogeneous field with values Hext = [Hx, Hy, Hz] for each cell and dimension [nx, ny, nz, 3]."""
        field = af.constant(0.0, nx, ny, nz, 3, dtype)
        field[:, :, :, 0]  = Hext[0]
        field[:, :, :, 1]  = Hext[1]
        field[:, :, :, 2]  = Hext[2]
        return field

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
    def __cinit__(self, unsigned nx, unsigned ny, unsigned nz, double dx, double dy, double dz):
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
        print (self._thisptr.nx)

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
    def dz(self):
        return self._thisptr.dz
    @dz.setter
    def dz(self, value):
        self._thisptr.dz=value



cdef class NonequiMesh:
    """
    Nonequispaced Mesh object.
    """
    def __init__(self, nx, ny, dx, dy, z_spacing):
        pass

    cdef cNonequiMesh* _thisptr
    cdef object owner # None if this is our own # From [1]
    #def __cinit__(self, int nx, int ny, double dx, double dy, vector[double] z_spacing):
    def __cinit__(self, unsigned nx, unsigned ny, double dx, double dy, z_spacing):
        cdef vector[double] z_spacing_cvec
        for val in z_spacing:
            z_spacing_cvec.push_back(val)
        self._thisptr = new cNonequiMesh(nx, ny, dx, dy, z_spacing)
        owner = None # see [1]
    cdef set_ptr(self, cNonequiMesh* ptr, owner):
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
    mean_m(i = None) : int
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
        if type(Ms) is af.Array:
            self._thisptr = new cState (deref(mesh._thisptr), <long int> addressof(Ms.arr), <long int> addressof(m.arr), <bool> verbose, <bool> mute_warning)
        else:
            self._thisptr = new cState (deref(mesh._thisptr), <double> Ms, <long int> addressof(m.arr), <bool> verbose, <bool> mute_warning)
        #af.device.lock_array(m_in)#This does not avoid memory corruption caused by double free
    def __dealloc__(self): # causes segfault on every cleanup
        del self._thisptr
        self._thisptr = NULL
    def __str__(self):
        """Convenience function, returns a string in the form 'state.t state.mean_mx() state.mean_my() state.mean_mz()'."""
        mx, my, mz = self.mean_m()
        return str(self.t) + " " + str(mx) + " " + str(my) + " " + str(mz)
        # return str(mx) + " " + str(my) + " " + str(mz)
        # return str(mx) + '\t' + str(my) + '\t' + str(mz)
        # return '<mx>=' + str(mx) + ', <my>=' + str(my) + ', <mz>=' + str(mz)

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
        Ms_field = array_from_addr(self._thisptr.wrapping_get_Ms_field())
        if Ms_field.is_empty():
            return self._thisptr.Ms
        else:
            return array_from_addr(self._thisptr.wrapping_get_Ms_field())
    @Ms.setter
    def Ms(self, value):
        if type(value) is af.Array:
            self._thisptr.set_Ms_field(<long int> addressof(value.arr))
        elif type(value) is float:
            self._thisptr.Ms = value
        else:
            raise RuntimeError("Ms.setter: value must be float or af.Array!")


    def write_vti(self, outputname):
        self._thisptr.write_vti( outputname.encode('utf-8'))
    def write_vti_atomistic(self, outputname):
        self._thisptr._vti_writer_atom( outputname.encode('utf-8'))
    def read_vti(self, inputname):
        self._thisptr._vti_reader( inputname.encode('utf-8'))

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
        return array_from_addr(self._thisptr.wrapping_get_Ms_field())
    @Ms_field.setter
    def Ms_field(self, Ms_field):
        self._thisptr.set_Ms_field(addressof(Ms_field.arr))

    @property
    def steps(self):
        return self._thisptr.steps
    def mean_M(self):
        """Spacial average <M>"""
        a = array_from_addr(self._thisptr.wrapping_mean_M_as_afarray())
        return a[0,0,0,0].scalar(), a[0,0,0,1].scalar(), a[0,0,0,2].scalar()

    def mean_m(self, i = None):
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
            return self._thisptr.mean_mx(), self._thisptr.mean_my(), self._thisptr.mean_mz()
        else:
            return self._thisptr.meani(i)
    def mean_mx(self):
        """Returns <mx>, i.e. average magnetisation in x-direction"""
        return self._thisptr.mean_mx()
    def mean_my(self):
        """Returns <my>, i.e. average magnetisation in y-direction"""
        return self._thisptr.mean_my()
    def mean_mz(self):
        """Returns <mz>, i.e. average magnetisation in z-direction"""
        return self._thisptr.mean_mz()


cdef class Stochastic_LLG:
    """
    Stochastic LLG Integrator.

    Parameters
    ----------
    alpha : float
        The unitless damping constant in the LLG equation
    T : float
        Temperature in Kelvin [K]
    dt : float
        Fixed time step in [s]
    terms : [FieldTerm]
        A python list constisting of FieldTerm objects s.a. ExchangeField or DemagField
    smode : str ("Heun")
        Integration scheme. Either "Heun" or "SemiHeun".
    """
    cdef cStochastic_LLG* _thisptr
    def __cinit__(self, alpha, T, dt, State state, terms=[], smode="Heun"):
        cdef vector[unique_ptr[cFieldterm]] vector_in
        if not terms:
            print("LLGIntegrator: no terms provided, please add some either by providing a list LLGIntegrator(terms=[...]) or calling add_terms(*args) after declaration.")
        else:
            for arg in terms:
                vector_in.push_back(unique_ptr[cFieldterm] (<cFieldterm*><size_t>arg._get_thisptr()))
            self._thisptr = new cStochastic_LLG (alpha, T, dt, deref(state._thisptr), move(vector_in), smode.encode('utf-8'))
    def step(self, State state):
        self._thisptr.step(deref(state._thisptr))
    def Eeff_in_J(self, State state):
        return self._thisptr.E(deref(state._thisptr))

cdef class LLGIntegrator:
    """
    LLGIntegrator(alpha, terms=[], mode="RKF45", hmin = 1e-15, hmax = 3.5e-10, atol = 1e-6, rtol = 1e-6)

    The LLGIntegrator object integrates a magnetization configuration according to the Landau–Lifshitz–Gilbert (LLG) equation

    Parameters
    ----------
    alpha : float
        The unitless damping constant in the LLG equation
    terms : [FieldTerm]
        A python list constisting of FieldTerm objects s.a. ExchangeField or DemagField
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
    Eeff_in_J(State) : float
        Calculates the micromagnetic energy of all terms for the magnetization state.m
    H_in_Apm(State) : af.array
        Returns the effective field H_in_Apm for the magnetization state.m
    add_terms(*args)
        Adds an FieldTerm object (s.a. ExchangeField) to be included in the effective field
    relax(State, precision, ncalcE, nprint)
        Relaxes the magnetization until the energy difference between ncalcE steps is less than precision


    Examples
    ----------
    llg = LLGIntegrator (alpha = 0.02, terms = terms, mode = "RKF45")
    llg.step(state)
    print(llg.Eeff_in_J(state))
    """

    cdef bool _alpha_is_array

    cdef cLLGIntegrator* _thisptr
    def __cinit__(self, alpha, terms, mode="RKF45", hmin = 1e-15, hmax = 3.5e-10, atol = 1e-6, rtol = 1e-6, dissipation_term_only = False):
        cdef vector[unique_ptr[cFieldterm]] vector_in
        if not terms:
            print("LLGIntegrator: no terms provided, please add some either by providing a list LLGIntegrator(terms=[...]) or calling add_terms(*args) after declaration.")
        else:
            for arg in terms:
                vector_in.push_back(unique_ptr[cFieldterm] (<cFieldterm*><size_t>arg._get_thisptr()))

            if hasattr(alpha, 'arr'):
                self._thisptr = new cLLGIntegrator (<long int> addressof(alpha.arr), move(vector_in), mode.encode('utf-8'), cController(hmin, hmax, atol, rtol), dissipation_term_only, False)
                self._alpha_is_array = True
            else:
                self._thisptr = new cLLGIntegrator (<double> alpha, move(vector_in), mode.encode('utf-8'), cController(hmin, hmax, atol, rtol), dissipation_term_only)
                self._alpha_is_array = False
    #def __dealloc__(self):
    #    # TODO maybe leads to segfault on cleanup, compiler warning eleminated by adding virtual destructor in adaptive_rk.hpp
    #    # NOTE is also problematic in minimizer class
    #    del self._thisptr
    #    self._thisptr = NULL
    def step(self, State state):
        self._thisptr.step(deref(state._thisptr))
    def Eeff_in_J(self, State state):
        return self._thisptr.E(deref(state._thisptr))
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr.h_addr(deref(state._thisptr)))
    def add_terms(self, *args):
        for arg in args:
            self._thisptr.llgterms.push_back(unique_ptr[cFieldterm] (<cFieldterm*><size_t>arg._get_thisptr()))
    def relax(self, State state, precision = 1e-10, ncalcE = 100, nprint = 1000, verbose = True):
        """
        relax(State state, precision = 1e-10, ncalcE = 100, nprint = 1000)
            Relaxes the magnetization until the energy difference between ncalcE steps is less than precision
        """
        self._thisptr.relax(deref(state._thisptr), precision, ncalcE, nprint, verbose)
    def integrate_dense(self, State state, double time_in_s, double write_every_dt_in_s, filename, verbose = False, append = False):
        """
            Integrate for time_in_s, using with dense output.
        """
        self._thisptr.integrate_dense(deref(state._thisptr), time_in_s, write_every_dt_in_s, filename.encode('utf-8'), verbose, append)
    @property
    def alpha(self):
        if self._alpha_is_array:
            return array_from_addr(self._thisptr.get_alpha_field_ptr())
        else:
            return self._thisptr.alpha
        # return self._thisptr.alpha
    @alpha.setter
    def alpha(self, value):
        if self._alpha_is_array:
            self._thisptr.set_alpha_field(addressof(value.arr))
        else:
            self._thisptr.alpha=value

    @property
    def accumulated_steps(self):
        return self._thisptr.accumulated_steps

        #cdef vector[unique_ptr[cFieldterm]] vector_in
        #for term in terms:
        #  vector_in.push_back(unique_ptr[cFieldterm] (<cFieldterm*><size_t>terms._get_thisptr()))
        #self._thisptr = new cLLGIntegrator (vector_in)

    #def print_stepsize(self):
    #  return self._thisptr.h_stepped_
    #def cpu_time(self):
    #  return self._thisptr.cpu_time()
    #def set_state0_alpha(self, value):
    #  self._thisptr.state0.material.alpha=value

cdef class StringMethod:
    """
    StringMethod method.
    """
    cdef cString* _thisptr
    #def __cinit__(self, State state, terms = [], n_interp=60, dt = 1e-13, LLGIntegrator llg):
    def __cinit__(self, State state, inputimages, n_interp, dt, LLGIntegrator llg):
        cdef vector[cState] vector_in
        if not inputimages:
            print("StringMethod: no States provided, please add at least two states for interpolation of initial path.")
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

            self._thisptr = new cString (deref(state._thisptr), move(vector_in), n_interp, dt, move(deref(llg._thisptr)))
            # Note: move(...llg._thisptr) prevents: "error: static assertion failed: result type must be constructible from value type of input range"
            # Note: we move vector_in also (is not strictly necessary here)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL

    def run(self, filepath, string_abort_rel_diff = 1e-12, string_abort_abs_diff = 1e-27, string_steps = 10000, every_string_to_vti = 50, verbose = True):
        return self._thisptr.run(filepath.encode('utf-8'), string_abort_rel_diff, string_abort_abs_diff, string_steps, every_string_to_vti, <bool> verbose)

cdef class FieldTerm:
    """ Field term base class, defines interface for calculating the effective field and micromagnetic energy."""
    def H_in_Apm(self, State state):
        """Effective field H in Ampere per meter [A/m]"""
        pass
    def H_in_T(self, State state):
        """Effective field H in Tesla [T]"""
        return Conversion.Apm_to_Tesla(self.H_in_Apm(state))
    def Energy_in_J(self, State state):
        """Micromagnetic energy in Joule [J]"""
        pass
    def Energy_in_eV(self, State state):
        """Micromagnetic energy in electronvolt [eV]"""
        return Conversion.J_to_eV(self.Energy_in_J(state))


cdef class DemagField(FieldTerm):
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
    def get_Nfft(self):
        return array_from_addr(self._thisptr.get_Nfft_ptr())
    ## Calculate energy contribution in [J]
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr

cdef class DemagFieldPBC(FieldTerm):
    """
    Demagnetization Field with Periodic Boundary Conditions.
    """
    cdef cDemagFieldPBC* _thisptr
    def __cinit__(self):
        self._thisptr = new cDemagFieldPBC ()
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class ExchangeField(FieldTerm):
    """
    Calculated field contribution of the exchange interaction, either by convolution (A is float) or as a sparse-matrix vector product (A is af.Array).
    The latter method properly handles jump conditions at material interfaces.

    Parameters
    ----------
    A : float or af.Array
        Exchange constant in  [J/m].
        If float: global exchange constant
        If af.Array: cell-wise defined exchange parameter.
    mesh : Mesh
        Discretization needed for setup of the sparse-matrix. Optional if A is float, must be provided otherwise.

    Examples
    ----------
    # (1) global value:
    A = 1e-15 # [J/m]
    SparseExchangeField(A)

    # (2) cell-wise value:
    mesh = Mesh(2, 1, 1, 1e-9, 1e-9, 1e-9)
    A_field = af.constant(0.0, mesh.nx, mesh.ny, mesh.nz, af.Dtype.f64)
    A_field[0] = 2e-15 # [J/m]
    A_field[1] = 3e-15 # [J/m]
    SparseExchangeField(A, mesh)
    """
    cdef bool _is_sparse
    cdef cExchangeField* _thisptr_conv
    cdef cSparseExchangeField* _thisptr_sparse
    def __cinit__(self, A, Mesh mesh = None, verbose = False):
        status_info = "\33[0;32mInfo:\33[0m"
        if mesh is None:
            if type(A) is af.Array:
                 raise RuntimeError("When passing A as af.Array, you must provide a mesh, i.e. ExchangeField(A, mesh).")
            elif type(A) is float:
                if verbose:
                    print(status_info, "ExchangeField.H_in_Apm using convolution")
                self._is_sparse = False
                self._thisptr_conv = new cExchangeField (<double> A)
            else:
                 raise RuntimeError("A must be float or af.Array!")
        elif type(mesh) is Mesh:
            self._is_sparse = True
            if verbose:
                print(status_info, "ExchangeField.H_in_Apm using sparse-matrix product")
            if type(A) is af.Array:
                self._thisptr_sparse = new cSparseExchangeField (<long int> addressof(A.arr), deref(mesh._thisptr), <bool> verbose)
            elif type(A) is float:
                # print("Note: you are using ExchangeField(float, Mesh) ")
                self._thisptr_sparse = new cSparseExchangeField (<double> A, deref(mesh._thisptr), <bool> verbose)
            else:
                raise RuntimeError("A must be either af.Array or float.")
        else:
            raise RuntimeError("mesh must be of type magnmaf.Mesh!")

    def __dealloc__(self):
        if self._is_sparse:
            del self._thisptr_sparse
            self._thisptr_sparse = NULL
        else:
            del self._thisptr_conv
            self._thisptr_conv = NULL

    def H_in_Apm(self, State state):
        if self._is_sparse:
            return array_from_addr(self._thisptr_sparse._pywrap_H_in_Apm(deref(state._thisptr)))
        else:
            return array_from_addr(self._thisptr_conv._pywrap_H_in_Apm(deref(state._thisptr)))

    def Energy_in_J(self, State state):
        if self._is_sparse:
            return self._thisptr_sparse.Energy_in_J(deref(state._thisptr))
        else:
            return self._thisptr_conv.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        if self._is_sparse:
            return self._thisptr_sparse.elapsed_eval_time()
        else:
            return self._thisptr_conv.elapsed_eval_time()
    def _get_thisptr(self):
        if self._is_sparse:
            return <size_t><void*>self._thisptr_sparse
        else:
            return <size_t><void*>self._thisptr_conv
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


cdef class SparseExchangeField(FieldTerm):
    cdef cSparseExchangeField* _thisptr
    def __cinit__(self, A, Mesh mesh, verbose = True):
        status_warning = "\33[1;31mWarning:\33[0m"
        print(status_warning, "SparseExchangeField is depricated and will be removed in the future, please use ExchangeField(A_field, mesh) instead.")
        if hasattr(A, 'arr'):
            self._thisptr = new cSparseExchangeField (<long int> addressof(A.arr), deref(mesh._thisptr), <bool> verbose)
        else:
            self._thisptr = new cSparseExchangeField (<double> A, deref(mesh._thisptr), <bool> verbose)
            # Note: use <bool_t> instead of <bool> in case of ambiguous overloading error: https://stackoverflow.com/questions/29171087/cython-overloading-no-suitable-method-found
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    def get_sparse_matrix(self):
        return array_from_addr(self._thisptr.pywrap_get_sparse_matrix_ptr())


cdef class ExchangeFieldPBC(FieldTerm):
    """
    Exchange Field with Periodic Boundary Condition (PBC). Implemented as sparse matrix.
    """
    cdef cExchangeFieldPBC* _thisptr
    def __cinit__(self, A, Mesh mesh, verbose = True):
        """
        Parameters
        ----------
        A : float
             Global exchange constant in [J/m]
        A : af.array
            cell-wise exchange constant defined at each cell, i.e. af.array of size [nx, ny, nz].
        mesh: mesh
        verbose: verbose switch
        """
        if hasattr(A, 'arr'):
            self._thisptr = new cExchangeFieldPBC (<long int> addressof(A.arr), deref(mesh._thisptr), <bool> verbose)
        else:
            self._thisptr = new cExchangeFieldPBC (<double> A, deref(mesh._thisptr), <bool> verbose)
            # Note: use <bool_t> instead of <bool> in case of ambiguous overloading error: https://stackoverflow.com/questions/29171087/cython-overloading-no-suitable-method-found
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr

cdef class BulkDMIExchangeField(FieldTerm):
    """
    Bulk Dzyaloshinskii–Moriya interaction (DMI) field.
    """
    cdef cBulkDMIExchangeField* _thisptr
    def __cinit__(self, D, A):
        """
        D : float
            Bulk DMI constant in units of [J/m]
        A : float
            Exchange constant in units of [J/m2]
        """
        self._thisptr = new cBulkDMIExchangeField (D, A)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()

cdef class DmiField(FieldTerm):
    """
    Interface Dzyaloshinskii–Moriya interaction (DMI) field.
    """
    cdef cDmiField* _thisptr
    def __cinit__(self, D, D_axis = [0., 0., 1.]):
        """
        Parameters
        ----------
        D : float
            DMI constant in units of [J/m^2]
        D : af.array (unsafe)
            cell-wise af.array of size [nx, ny, nz]. Caution: does not yet support jump conditions
        D_axis : [float, float, float]
            DMI orientation axis
            defaults to [0., 0., 1.]
        """
        if hasattr(D, 'arr'):
            self._thisptr = new cDmiField (<long int> addressof(D.arr), <double> D_axis[0], <double> D_axis[1], <double> D_axis[2])
        else:
            self._thisptr = new cDmiField (<double> D, <double> D_axis[0], <double> D_axis[1], <double> D_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr

cdef class DMI_D2d_Field(FieldTerm):
    """
    Dzyaloshinskii–Moriya interaction (DMI) for materials of crystallographic class D2d.
    """
    cdef cDMI_D2d_Field* _thisptr
    def __cinit__(self, D, PBC = False):
        """
        Parameters
        ----------
        D : float
            DMI constant in units of [J/m^2]
        PBC : bool
            if true, use periodic boundary conditions.
        """
        self._thisptr = new cDMI_D2d_Field (<double> D, <bool> PBC)
    @property
    def D_in_J_per_m2_(self):
        return self._thisptr.D_in_J_per_m2_
    @D_in_J_per_m2_.setter
    def D_in_J_per_m2_(self, value):
      self._thisptr.D_in_J_per_m2_=value

    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr

cdef class NonequiDemagField(FieldTerm):
    cdef cNonequiDemagField* _thisptr
    def __cinit__(self, NonequiMesh mesh, verbose = True, caching = True, nthreads : int = 8):
        self._thisptr = new cNonequiDemagField (deref(mesh._thisptr), <bool> verbose, <bool> caching, nthreads)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def get_Nfft(self):
        return array_from_addr(self._thisptr.get_Nfft_ptr())
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
        return <size_t><void*>self._thisptr


cdef class NonequiExternalField(FieldTerm):
    cdef cNonequiExternalField* _thisptr
    def __cinit__(self, NonequiMesh mesh, array_in):
        self._thisptr = new cNonequiExternalField (deref(mesh._thisptr), addressof(array_in.arr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    # def set_homogeneous_field(self, x, y, z):
    #         self._thisptr.set_homogeneous_field(x, y, z)
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class NonequiExchangeField(FieldTerm):
    cdef cNonequiExchangeField* _thisptr
    def __cinit__(self, NonequiMesh mesh, A, verbose = True):
        if hasattr(A, 'arr'):
            self._thisptr = new cNonequiExchangeField (deref(mesh._thisptr), <long int> addressof(A.arr), <bool> verbose)
        else:
            self._thisptr = new cNonequiExchangeField (deref(mesh._thisptr), <double> A, <bool> verbose)
            # Note: use <bool_t> instead of <bool> in case of ambiguous overloading error: https://stackoverflow.com/questions/29171087/cython-overloading-no-suitable-method-found
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
        return <size_t><void*>self._thisptr

cdef class CubicAnisotropyField(FieldTerm):
    """
    Cubic magneto-crystalline anisotropy.
    """
    cdef cCubicAnisotropyField* _thisptr
    def __cinit__(self, Kc1, Kc2 = 0, Kc3 = 0, c1 = [1, 0, 0], c2 = [0, 1, 0]):
        """
        Parameters
        ----------
        Kc1 : float or array [nx, ny, nz, 1]
            1st order cubic anisotropy constant in [J/m^3]
        Kc2 : float or array [nx, ny, nz, 1]
            2nd order cubic anisotropy constant in [J/m^3]
        Kc3 : float or array [nx, ny, nz, 1]
            3rd order cubic anisotropy constant in [J/m^3]
        c1 : [float, float, float] or array [nx, ny, nz, 3]
            unit vector indicating anisotropy direction
            defaults to [1., 0., 0.]
        c2 : [float, float, float] or array [nx, ny, nz, 3]
            unit vector indicating anisotropy direction, must be choosen orthogonal to c1.
            defaults to [0., 1., 0.]
        """
        if hasattr(Kc1, 'arr') and hasattr(Kc2, 'arr') and hasattr(Kc3, 'arr') and hasattr(c1, 'arr') and hasattr(c2, 'arr'):
            self._thisptr = new cCubicAnisotropyField (<long int> addressof(Kc1.arr), <long int> addressof(Kc2.arr), <long int> addressof(Kc3.arr), <long int> addressof(c1.arr), <long int> addressof(c2.arr))
        else:
            self._thisptr = new cCubicAnisotropyField (<double> Kc1, <double> Kc2, <double> Kc3, <double> c1[0], <double> c1[1], <double> c1[2], <double> c2[0], <double> c2[1], <double> c2[2])
    # tried this i.o.t. wrap std::array<double, 3>:
    #def __cinit__(self, Kc1, Kc2, Kc3, double [:] c1, double [:] c2):
    #    cdef carray_d3 *c1arr = <carray_d3 *>(&c1[0])
    #    cdef carray_d3 *c2arr = <carray_d3 *>(&c2[0])
    #    self._thisptr = new cCubicAnisotropyField (Kc1, Kc2, Kc3, c1arr[0], c2arr[0])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
        return <size_t><void*>self._thisptr


cdef class UniaxialAnisotropyField(FieldTerm):
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
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    @property
    def Ku1(self):
        return self._thisptr.Ku1
    @Ku1.setter
    def Ku1(self, value):
      self._thisptr.Ku1=value

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


cdef class NonequiUniaxialAnisotropyField(FieldTerm):
    cdef cNonequiUniaxialAnisotropyField* _thisptr
    def __cinit__(self, NonequiMesh mesh, Ku1, Ku1_axis = [0, 0, 1]):
        if hasattr(Ku1, 'arr') and hasattr(Ku1_axis, 'arr'):
            self._thisptr = new cNonequiUniaxialAnisotropyField (deref(mesh._thisptr), <long int> addressof(Ku1.arr), <long int> addressof(Ku1_axis.arr))
        elif hasattr(Ku1, 'arr') :
            self._thisptr = new cNonequiUniaxialAnisotropyField (deref(mesh._thisptr), <long int> addressof(Ku1.arr), <double> Ku1_axis[0], <double> Ku1_axis[1], <double> Ku1_axis[2])
        elif hasattr(Ku1_axis, 'arr'):
            self._thisptr = new cNonequiUniaxialAnisotropyField (deref(mesh._thisptr), <double> Ku1, <long int> addressof(Ku1_axis.arr))
        else:
            self._thisptr = new cNonequiUniaxialAnisotropyField (deref(mesh._thisptr), <double> Ku1, <double> Ku1_axis[0], <double> Ku1_axis[1], <double> Ku1_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
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


cdef class AtomisticDipoleDipoleField(FieldTerm):
    cdef cAtomisticDipoleDipoleField* _thisptr
    def __cinit__(self, Mesh mesh):
        self._thisptr = new cAtomisticDipoleDipoleField (deref(mesh._thisptr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticUniaxialAnisotropyField(FieldTerm):
    cdef cAtomisticUniaxialAnisotropyField* _thisptr
    def __cinit__(self, double K_atom, K_atom_axis = [0., 0., 1.]):
    #def __cinit__(self, double K_atom, double K_atom_axis_x, double K_atom_axis_y, double K_atom_axis_z):
        self._thisptr = new cAtomisticUniaxialAnisotropyField (K_atom, K_atom_axis[0], K_atom_axis[1], K_atom_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticExchangeField(FieldTerm):
    cdef cAtomisticExchangeField* _thisptr
    def __cinit__(self, double J_atom):
        self._thisptr = new cAtomisticExchangeField (J_atom)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticDmiField(FieldTerm):
    cdef cAtomisticDmiField* _thisptr
    def __cinit__(self, double D_atom, D_atom_axis = [0., 0., -1.]):
        self._thisptr = new cAtomisticDmiField (D_atom, D_atom_axis[0], D_atom_axis[1], D_atom_axis[2])
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class ExternalField(FieldTerm):
    cdef cExternalField* _thisptr
    def __cinit__(self, array_in):
        self._thisptr = new cExternalField (addressof(array_in.arr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def set_homogeneous_field(self, x, y, z):
            self._thisptr.set_homogeneous_field(x, y, z)
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class AtomisticExternalField(FieldTerm):
    cdef cAtomisticExternalField* _thisptr
    def __cinit__(self, array_in):
        self._thisptr = new cAtomisticExternalField (addressof(array_in.arr))
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def cpu_time(self):
        return self._thisptr.elapsed_eval_time()
    def set_homogeneous_field(self, x, y, z):
            self._thisptr.set_homogeneous_field(x, y, z)
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class SpinTransferTorqueField(FieldTerm):
    """
    Slonczewski Spin Transfer Torque model

    Parameters
    ----------
    pol : af.array [nx, ny, nz, 3]
        Unit length polarization magnetization vector at each cell in [a.u.].
    eta_damp : float
        Damping like constant in [a.u.]
    eta_field : float
        Field like constant in [a.u.]
    j_e : float
        Current in [A]
    fl_thickness : float
        Free layer thickness in [m]
    """
    cdef cSpinTransferTorqueField* _thisptr
    def __cinit__(self, pol, eta_damp : float,  eta_field : float, j_e : float, fl_thickness : float):
        self._thisptr = new cSpinTransferTorqueField (addressof(pol.arr), eta_damp, eta_field, j_e, fl_thickness)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr
    @property
    def eta_damping(self):
        return self._thisptr.eta_damping
    @eta_damping.setter
    def eta_damping(self, value):
        self._thisptr.eta_damping = value

    @property
    def eta_field(self):
        return self._thisptr.eta_field
    @eta_field.setter
    def eta_field(self, value):
        self._thisptr.eta_field = value

    @property
    def j_e(self):
        return self._thisptr.j_e
    @j_e.setter
    def j_e(self, value):
        self._thisptr.j_e = value

    @property
    def fl_thickness(self):
        return self._thisptr.fl_thickness
    @fl_thickness.setter
    def fl_thickness(self, value):
        self._thisptr.fl_thickness = value


    # @property
    # def polarization_field(self):
    #     return array_from_addr(self._thisptr.polarization_field.get_array_addr())
    # @polarization_field.setter
    # def polarization_field(self, array):
    #     self._thisptr.polarization_field.set_array(addressof(array.arr))


cdef class RKKYExchangeField(FieldTerm):
    cdef cRKKYExchangeField* _thisptr
    def __cinit__(self, rkky_values, exchange_values, Mesh mesh, rkky_indices = None, verbose = True):
        if rkky_indices is None: # setting default zeros as passing empty array is not possible
            rkky_indices = af.constant(0, mesh.nx, mesh.ny, mesh.nz, 3, dtype=af.Dtype.u32)
        self._thisptr = new cRKKYExchangeField (addressof(rkky_values.arr), addressof(exchange_values.arr), deref(mesh._thisptr), addressof(rkky_indices.arr), <bool> verbose)
    def __dealloc__(self):
        del self._thisptr
        self._thisptr = NULL
    def H_in_Apm(self, State state):
        return array_from_addr(self._thisptr._pywrap_H_in_Apm(deref(state._thisptr)))
    def Energy_in_J(self, State state):
        return self._thisptr.Energy_in_J(deref(state._thisptr))
    def _get_thisptr(self):
            return <size_t><void*>self._thisptr


cdef class LBFGS_Minimizer:
    """
    The LBFGS_Minimizer object implements an Limited-memory Broyden–Fletcher–Goldfarb–Shanno (LBFGS) algorithm, minimizing a magentization configuratio with respect to the micromagnetic energy. This energy is obtained via the effective field given by the FieldTerm objects.

    Parameters
    ----------
    terms : [FieldTerm]
        List of FieldTerm objects
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
    add_terms([FieldTerm])
        Adds a list of FieldTerm objects to the existing terms
    delete_last_term()
        Delets the lst FieldTerm object in the list

    Examples
    --------
    minimizer = Minimizer(terms, tol = 1e-6, maxiter = 230)
    minimier.minimize(state)
    """
    cdef cLBFGS_Minimizer* _thisptr
    def __cinit__(self, terms=[], tol = 1e-6, maxiter = 230, verbose = 0):
        cdef vector[unique_ptr[cFieldterm]] vector_in
        if not terms:
            print("cLBFGS_Minimizer: no terms provided, please add some either by providing a list terms=[...] or calling add_terms(*args)")
        else:
            for arg in terms:
                vector_in.push_back(unique_ptr[cFieldterm] (<cFieldterm*><size_t>arg._get_thisptr()))
            self._thisptr = new cLBFGS_Minimizer (move(vector_in), tol, maxiter, verbose)
            # Note: std::cout is not handled and leads to segfault on some installs, hence c++ backend uses prinft

    #TODO
    #  def __dealloc__(self): #causes segfault on GTO in cleanup
    #      del self._thisptr
    #      self._thisptr = NULL
    def add_terms(self, *args):
        for arg in args:
            self._thisptr.fieldterms_.push_back(unique_ptr[cFieldterm] (<cFieldterm*><size_t>arg._get_thisptr()))
    def delete_last_term(self):
        self._thisptr.fieldterms_.pop_back()
    def minimize(self, State state):
        return self._thisptr.Minimize(deref(state._thisptr))


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
    e_neg : float
        [C] negative signed elementary charge e
    e_abs : float
        [C] absolute value of elementary charge e
    kb : float
        [J/K] Boltzmann constant kb
    hbar : float
        [J s] reduced Planck constant
    """
    mu0 = 4e-7 * pi

    gamma = 1.760859644e11 * mu0

    mu_b = 9.274009994e-24

    e_neg = - 1.6021766208e-19

    e_abs = 1.6021766208e-19

    kb = 1.38064852e-23

    hbar = 1.0545718e-34
