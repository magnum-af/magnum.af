cmake_minimum_required(VERSION 3.5)
project(magnum.af
    VERSION 1.0.0
    DESCRIPTION "A finite differences GPU-accelerated micromagnetic and atomistic simulation software" )

message(STATUS "Setting up project ${PROJECT_NAME} version ${PROJECT_VERSION}.")

# options consumed in ${magnumaf_core_dir_name}:
option(EXAMPLES "Build examples." ON)
option(TESTS "Build unit tests." ON)
option(DEBUG "Build with debugging flag '-g'." OFF)
option(SETINSTALL "Building shared library, configuring install directories." ON)
option(PACK_AS_DEB "Pack as .deb file instead of package script and .tar files." ON)
option(CUDA "Build CUDA backend if available." ON)
option(OPENCL "Build OPENCL backend if available." ON)
option(DOCU "Build doxygen documentation if installed." ON)

# option consumed in python/
option(PYTHON "Build python bindings" ON)

if(TESTS)
  message(STATUS "Enable unit testing (option 'TESTING').")
  enable_testing() # enables CTest with 'make test' and 'ctest' commands
  # note: enable must come before we include tests below
endif()

set(magnumaf_core_dir_name "magnumaf") # note: variable is used in python/CMakeLists.txt
add_subdirectory(${magnumaf_core_dir_name})

# build python bindings
if(PYTHON)
  message(STATUS "Building python bindings (option 'PYTHON').")
  add_subdirectory(python)
endif()

if(DOCU)
  add_subdirectory(docu)
endif()
