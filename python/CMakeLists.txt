CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(magnum.af)# this adapts compiler flags from parent

# building with python3 only
set(PYTHON_EXECUTABLE "/usr/bin/python3") # equivalent to cmake -DPYTHON_EXECUTABLE=/usr/bin/python3 # hint from https://github.com/pytorch/ELF/issues/98
#Note: if python3 is not aviable use:# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing") # in response to "dereferencing type-punned pointer will break strict-aliasing rules" warning # hint: https://github.com/cython/cython/issues/958

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)
# Include the CMake script UseCython.cmake.  This defines add_cython_module().
# Instruction for use can be found at the top of CMakeModules/UseCython.cmake.
include( UseCython )
##Currently not needed, but maybe interesting at some point:
# With CMake, a clean separation can be made between the source tree and the
# build tree.  When all source is compiled, as with pure C/C++, the source is
# no-longer needed in the build tree.  However, with pure *.py source, the
# source is processed directly.  To handle this, we reproduce the availability
# of the source files in the build tree.

#add_custom_target( ReplicatePythonSourceTree ALL ${CMAKE_COMMAND} -P
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake
#  ${CMAKE_CURRENT_BINARY_DIR}
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
#include_directories( ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/include ) # Warning: commenting this out leads to "fatal error: ../../src/LLGTerm.hpp: No such file or directory compilation terminated."
# If the pyx file is a C++ file, we should specify that here.
set_source_files_properties(
  ${magnum.af_SOURCE_DIR}/magnumaf.pyx
  #Additional files would come here as follows:
  #${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/FurtherFiles.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )

FILE(GLOB SOURCE_FOR_PYBINDING "../src/*.cpp" "../src/*/*.cpp")
cython_add_module_w_AF_VTK( magnumaf magnumaf.pyx ${SOURCE_FOR_PYBINDING})
# Further modules would come with cython_add_module_w_AF_VTK(target target.pyx [sources])

## Linking with magnum.af_CORE_SHARED would be useful for avoiding second source build, but cython is currently now working with shared library
#FIND_PACKAGE(ArrayFire REQUIRED)
#INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS})# Include the ArrayFire hreaders
#ADD_LIBRARY(magnum.af_CORE_OBJECT OBJECT ${SOURCE_FOR_PYBINDING}) # TODO error with python when not OBJECT
#cython_add_module_w_AF_VTK( magnumaf magnumaf.pyx $<TARGET_OBJECTS:magnum.af_CORE_OBJECT>)