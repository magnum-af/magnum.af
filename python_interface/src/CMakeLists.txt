# If the pyx file is a C++ file, we should specify that here.
FIND_PACKAGE(ArrayFire REQUIRED)
#If VTK_DIR not set, you may set: set(VTK_DIR "/usr/include/vtk-6.2/")
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
#FIND_PACKAGE(OpenCL)
#IF(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
#    MESSAGE(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL benchmark")
#ENDIF()
set_source_files_properties(
#  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/rect.pyx
#  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/pi_calculator.pyx
  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/interface.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )


# Multi-file cython modules do not appear to be working at the moment.
#cython_add_module( rect rect.pyx Rectangle.cxx )

#cython_add_module( pi_calculator pi_calculator.pyx PiCalculator.cxx )

cython_add_module_w_AF_VTK( interface interface.pyx ../../src/mesh.cpp ../../src/state.cpp ../../src/micro_demag.cpp ../../src/func.cpp ../../src/llg.cpp ../../src/micro_exch.cpp ../../src/zee.cpp ../../src/vtk_writer.cpp)
#set( PYTHON_LIBRARIES ${PYTHON_LIBRARIES} ${ArrayFire_OpenCL_LIBRARIES} ${VTK_LIBRARIES})
