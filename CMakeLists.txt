cmake_minimum_required(VERSION 3.2)

project(magnum.af VERSION 0.1.0)
message(STATUS "Setting up project ${PROJECT_NAME} version ${PROJECT_VERSION}.")

option(PYTHON "Build python bindings" ON)
# consumed in ${magnumaf_core_dir_name}:
option(EXAMPLES "Build examples." ON)
option(TESTS "Build unit tests." ON)
option(DEBUG "Build with debugging flag '-g'." OFF)
option(SETINSTALL "Building shared library, configuring install directories." ON)
option(PACK_AS_DEB "Pack as .deb file instead of package script and .tar files." ON)
option(CUDA "Build CUDA backend." ON)
option(OPENCL "Build OPENCL backend." ON)

set(magnumaf_core_dir_name "magnumaf") # note: variable is used in python/CMakeLists.txt
add_subdirectory(${magnumaf_core_dir_name})

if(TESTS)
  message(STATUS "Enable unit testing (option 'TESTING').")
  enable_testing() # enables CTest with 'make test' and 'ctest' commands
endif()

# build python bindings
if(PYTHON)
  message(STATUS "Building python bindings (option 'PYTHON').")
  add_subdirectory(python)
endif()
