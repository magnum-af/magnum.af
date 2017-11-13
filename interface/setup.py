from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'PTH',
  ext_modules=[
    Extension("pth_mag",
              # Note, you can link against a c++ library
              # instead of including the source
              libraries = ["afcuda"],
              #TODO libraries = ["afcpu"],
              sources=["interface.pyx", "../src/mesh.cpp","../src/state.cpp"],
              language="c++"),
    ],
  cmdclass = {'build_ext': build_ext},
)
