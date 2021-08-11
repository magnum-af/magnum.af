project(${PROJECT_NAME} VERSION ${PROJECT_VERSION}) # NOTE: setting project
                                                    # again prevents 'operation
                                                    # only allowed in c++ file
                                                    # error"
# sets magnum.af_VERSION_XXX and PROJECT_VERSION_XXX variables, with XXX =
# MAJOR,MINOR,PATCH,TWEAK

# This is redundant, as mentioned in c++ part already:
include_directories(
  ${CMAKE_SOURCE_DIR}/${magnumaf_core_dir_name}/include/magnumaf
  ${CMAKE_SOURCE_DIR}/${magnumaf_core_dir_name}/include/magnumaf/field_terms
  ${CMAKE_SOURCE_DIR}/${magnumaf_core_dir_name}/include/magnumaf/solvers
  ${CMAKE_SOURCE_DIR}/${magnumaf_core_dir_name}/include/magnumaf/util
  ${CMAKE_SOURCE_DIR}/${magnumaf_core_dir_name}/include/magnumaf/integrators
  ${CMAKE_BINARY_DIR}/${magnumaf_core_dir_name}/include/magnumaf/ # for
                                                                  # configured
                                                                  # headers
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
if(DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()
# disables warning caused by cython compiled file
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-conversion")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-cast-function-type")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wimplicit-fallthrough=0")

# building with python3 only
set(PYTHON_EXECUTABLE "/usr/bin/python3")
# equivalent to cmake -DPYTHON_EXECUTABLE=/usr/bin/python3 # hint from
# https://github.com/pytorch/ELF/issues/98 Note: if python3 is not aviable use:#
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing") # in response
# to "dereferencing type-punned pointer will break strict-aliasing rules"
# warning # hint: https://github.com/cython/cython/issues/958

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)
# Include the CMake script UseCython.cmake.  This defines add_cython_module().
# Instruction for use can be found at the top of CMakeModules/UseCython.cmake.
include(UseCython)
# Currently not needed, but maybe interesting at some point: With CMake, a clean
# separation can be made between the source tree and the build tree.  When all
# source is compiled, as with pure C/C++, the source is no-longer needed in the
# build tree.  However, with pure *.py source, the source is processed directly.
# To handle this, we reproduce the availability of the source files in the build
# tree.

# add_custom_target( ReplicatePythonSourceTree ALL ${CMAKE_COMMAND} -P
# ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake
# ${CMAKE_CURRENT_BINARY_DIR} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
# include_directories( ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/include ) # Warning:
# commenting this out leads to "fatal error: ../../src/field_term.hpp: No such
# file or directory compilation terminated." If the pyx file is a C++ file, we
# should specify that here.
set_source_files_properties(
  ${magnum.af_SOURCE_DIR}/magnumaf.pyx
  # Additional files would come here as follows:
  # ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/FurtherFiles.pyx
  PROPERTIES CYTHON_IS_CXX TRUE)

# FILE(GLOB SOURCE_FOR_PYBINDING "../src/*.cpp" "../src/*/*.cpp")
# cython_add_module_w_AF_VTK( magnumaf magnumaf.pyx ${SOURCE_FOR_PYBINDING})
cython_add_module_w_af_vtk(magnumaf magnumaf.pyx
                           $<TARGET_OBJECTS:magnumaf_object_lib>)

if(SETINSTALL)
  # 'make install' to the correct locations (provided by GNUInstallDirs).
  install(
    TARGETS magnumaf
    EXPORT magnumafConfig
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # This is for Windows
  )

  set_target_properties(magnumaf PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE
  )# prevent seting runtime path RPATH to ""
  # from:
  # https://stackoverflow.com/questions/32469953/why-is-cmake-designed-so-that-it-removes-runtime-path-when-installing
endif()

# MESSAGE( STATUS "------> CMAKE_CXX_FLAGS= ${CMAKE_CXX_FLAGS}" )

# Further modules would come with cython_add_module_w_AF_VTK(target target.pyx
# [sources])

# Linking with magnum.af_CORE_SHARED would be useful for avoiding second source
# build, but cython is currently now working with shared library
# FIND_PACKAGE(ArrayFire REQUIRED)
# INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS})# Include the ArrayFire hreaders
# ADD_LIBRARY(magnumaf_object_lib OBJECT ${SOURCE_FOR_PYBINDING}) # TODO error
# with python when not OBJECT cython_add_module_w_AF_VTK( magnumaf magnumaf.pyx
# $<TARGET_OBJECTS:magnumaf_object_lib>)
