CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(magnum.af)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wuseless-cast")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast")

# Find the ArrayFire package.
FIND_PACKAGE(ArrayFire REQUIRED)
# If ArrayFire is found, the following variables will be defined:
#
# ArrayFire_INCLUDE_DIRS    - Location of ArrayFire's include directory.
# ArrayFire_LIBRARIES       - Location of ArrayFire's libraries. This will default
#                             to a GPU backend if one is found.
# ArrayFire_FOUND           - True if ArrayFire has been located
#
# You may provide a hint to where ArrayFire's root directory may be located
# by setting ArrayFire_DIR.
#
# ----------------------------------------------------------------------------
#
# ArrayFire_CPU_FOUND        - True of the ArrayFire CPU library has been found.
# ArrayFire_CPU_LIBRARIES    - Location of ArrayFire's CPU library, if found
# ArrayFire_CUDA_FOUND       - True of the ArrayFire CUDA library has been found.
# ArrayFire_CUDA_LIBRARIES   - Location of ArrayFire's CUDA library, if found
# ArrayFire_OpenCL_FOUND     - True of the ArrayFire OpenCL library has been found.
# ArrayFire_OpenCL_LIBRARIES - Location of ArrayFire's OpenCL library, if found

#VKT
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})


# Include the ArrayFire hreaders
INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS})

# Get a list of all of the source files in the directory:
FILE(GLOB SOURCE "*.cpp" "*/*.cpp")

# There are several ways of compiling source code in CMake. In most cases you
# specify the source files to an ADD_EXCUTABLE call. Because we intend to
# link against multiple ArrayFire backends, this approach would result in
# every source file being re-compiled multiple times.
# Here we take a different approach to compilation. We build all of the source
# files into a CMake object library and then link this library against the
# different ArrayFire backends.

ADD_LIBRARY(magnumafcpp SHARED ${SOURCE})
ADD_LIBRARY(magnum.af_CORE_OBJECT OBJECT ${SOURCE})

# following: https://rix0r.nl/blog/2015/08/13/cmake-guide/

include_directories( ../include ../include/llg_terms ../include/solvers ../include/util ../include/integrators )

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
#target_include_directories(magnumafcpp PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include> # TODO remove /..
#    $<INSTALL_INTERFACE:../include>
#    PRIVATE ../src) # TODO remove /..
#


# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS magnumafcpp EXPORT magnumafcppConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
#TODO with correct headers # install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
install(EXPORT magnumafcppConfig DESTINATION share/magnumaf/cmake)

# This makes the project importable from the build directory
export(TARGETS magnumafcpp FILE magnumafcppConfig.cmake)
