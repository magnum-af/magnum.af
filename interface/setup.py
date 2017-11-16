from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'PTH',
  ext_modules=[
    Extension("pth_mag",
              # Note, you can link against a c++ library
              # instead of including the source
              #libraries = ["afcuda"],
              libraries = ["afcpu"],
              #sources=["interface.pyx", "../src/mesh.cpp","../src/state.cpp","../src/micro_demag.cpp","../src/func.cpp"],
              sources=["interface.pyx", "../src/mesh.cpp","../src/state.cpp","../src/micro_demag.cpp","../src/func.cpp","../src/llg.cpp","../src/micro_exch.cpp","../src/zee.cpp"],
              extra_compile_args=['-std=gnu++14','-O3'],
              language="c++"),
    ],
  cmdclass = {'build_ext': build_ext},
)
