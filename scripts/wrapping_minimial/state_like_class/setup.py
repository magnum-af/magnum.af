from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'Demos',
  ext_modules=[
    Extension("test",
              libraries = ["afcuda"],
              library_dirs = ["/usr/local/arrayfire/lib"],
              include_dirs = ["/usr/local/arrayfire/include"],
              sources=["test.pyx", "cpp_test.cpp"],
              language="c++"),
    ],
  cmdclass = {'build_ext': build_ext},
)

