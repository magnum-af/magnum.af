CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(pth-mag)
#Note in new projects, set VTK_DIR by export VTK_DIR=/home/.../VTK-build:$VTK_DIR

# Compile with C++14 support.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
    message(STATUS "Error: The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler or update the existing one.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# indicate we wish for the binaries to be written to the bin directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
# we have custom CMake configuration files, tell CMake to use them:
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
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


# build items in subdirectories
add_subdirectory(src)
