from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'Wrap',
  ext_modules=[
    Extension("wrap",
              #libraries = ["afcuda"],
              libraries = ["afopencl"],
              #libraries = ["afcpu"],
              sources=["interface.pyx", "source.cpp"],
              extra_compile_args=['-std=gnu++14', '-O3', '-I/opt/arrayfire/include'],
              extra_link_args=['-L/opt/arrayfire/lib'],
              language="c++"),
    ],
  cmdclass = {'build_ext': build_ext},
)
