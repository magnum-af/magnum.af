from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

#check: https://github.com/r9y9/pypcl/blob/master/setup.py for vtk settings

extensions = [
    Extension("magnum_af_cython_setup", ["magnumaf.pyx"],
        include_dirs = ["/opt/include", "/usr/include/vtk-7.1"],
        #libraries = ["afcuda"],
        #libraries = ["afopencl"],
        libraries = ["afcpu"],
        library_dirs = ["/opt/lib", "/usr/include/vtk-7.1"], 
        extra_compile_args=["-std=c++11"],
        )
]
setup(
    name="magnum_af_cython_setup",
    ext_modules=cythonize(extensions),
)
