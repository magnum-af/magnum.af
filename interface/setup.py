from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


setup(
  name = 'PTH',
  ext_modules=[
    Extension("pth_mag",
              # Note, you can link against a c++ library
              # instead of including the source
              #libraries = ["afcuda"],
              #libraries = ["afcpu"],
              #libraries = ["afopencl"],
              #libraries = ["afopencl","vtkWrappingTools-6.2"],
              libraries = ["afopencl"],
              #libraries = ["afopencl","vtkCommonCore-6.2"],
              #sources=["interface.pyx", "../src/mesh.cpp","../src/state.cpp","../src/micro_demag.cpp","../src/func.cpp"],
              #sources=["interface.pyx", cythonize("src/*.cpp", include_path=[..])], # include_path = [../src/]],
              #sources=["interface.pyx",cythonize("*.cpp", include_path = ["../src"])], # include_path = [../src/]],
              #sources=["interface.pyx", "../src/mesh.cpp","../src/state.cpp","../src/micro_demag.cpp","../src/func.cpp","../src/llg.cpp","../src/micro_exch.cpp","../src/zee.cpp"],
              sources=["interface.pyx", "../src/mesh.cpp","../src/state.cpp","../src/micro_demag.cpp","../src/func.cpp","../src/llg.cpp","../src/micro_exch.cpp","../src/zee.cpp"],
              #sources=["interface.pyx", "../src/mesh.cpp","../src/state.cpp","../src/micro_demag.cpp","../src/func.cpp","../src/llg.cpp","../src/micro_exch.cpp","../src/zee.cpp","../src/vtk_writer.cpp"],
              #extra_compile_args=['-std=gnu++14','-O3'],
              #extra_compile_args=['-std=gnu++14','-O3'],
              extra_compile_args=['-std=gnu++14','-O3'],
              #extra_compile_args=['-std=gnu++14','-O3','-I/usr/include/vtk-6.2'],
              #extra_compile_args=['-std=gnu++14','-O3','-DVTK_DIR:PATH=/home/pth/programs/VTK_build','-I/usr/include/vtk-6.2'],
              #extra_compile_args=['-std=gnu++14','-O3','-DVTK_DIR:PATH=/home/pth/programs/VTK_build','-I/usr/local/include/vtk-8.1'],
              #extra_compile_args=['-std=gnu++14','-O3','-DVTK_DIR:PATH=/home/pth/programs/VTK_build'],
              language="c++"),
    ],
  cmdclass = {'build_ext': build_ext},
)
