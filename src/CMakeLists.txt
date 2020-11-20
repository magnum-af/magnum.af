set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wuseless-cast")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast")

# Find the ArrayFire package.
find_package(ArrayFire REQUIRED)
# If ArrayFire is found, the following variables will be defined:
#
# ArrayFire_INCLUDE_DIRS    - Location of ArrayFire's include directory.
# ArrayFire_LIBRARIES       - Location of ArrayFire's libraries. This will
# default to a GPU backend if one is found. ArrayFire_FOUND           - True if
# ArrayFire has been located
#
# You may provide a hint to where ArrayFire's root directory may be located by
# setting ArrayFire_DIR.
#
# ----------------------------------------------------------------------------
#
# ArrayFire_CPU_FOUND        - True of the ArrayFire CPU library has been found.
# ArrayFire_CPU_LIBRARIES    - Location of ArrayFire's CPU library, if found
# ArrayFire_CUDA_FOUND       - True of the ArrayFire CUDA library has been
# found. ArrayFire_CUDA_LIBRARIES   - Location of ArrayFire's CUDA library, if
# found ArrayFire_OpenCL_FOUND     - True of the ArrayFire OpenCL library has
# been found. ArrayFire_OpenCL_LIBRARIES - Location of ArrayFire's OpenCL
# library, if found

# VKT
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Include the ArrayFire headers (SYSTEM option disables warnings from afheaders)
include_directories(SYSTEM ${ArrayFire_INCLUDE_DIRS})

# Get a list of all of the source files in the directory:
file(GLOB SOURCE "*.cpp" "*/*.cpp" "*/*/*.cpp")

# There are several ways of compiling source code in CMake. In most cases you
# specify the source files to an ADD_EXCUTABLE call. Because we intend to link
# against multiple ArrayFire backends, this approach would result in every
# source file being re-compiled multiple times. Here we take a different
# approach to compilation. We build all of the source files into a CMake object
# library and then link this library against the different ArrayFire backends.

# adding OBJECT library instead of SHARED such that the python binding doesn't
# have to recompile all sources:
add_library(magnumaf_object_lib OBJECT ${SOURCE})
# explicitly turning PIC on (in case the compiler flag is not set):
set_property(TARGET magnumaf_object_lib PROPERTY POSITION_INDEPENDENT_CODE ON)

if(SETINSTALL)
  # installing magnumafcpp lib, adapted from:
  # https://rix0r.nl/blog/2015/08/13/cmake-guide/
  add_library(magnumafcpp SHARED $<TARGET_OBJECTS:magnumaf_object_lib>)

  # Define headers for this library. PUBLIC headers are used for compiling the
  # library, and will be added to consumers' build paths.
  target_include_directories(
    magnumafcpp
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include> # TODO remove
                                                                # /..
      $<INSTALL_INTERFACE:include>
    PRIVATE ../src) # TODO remove /..

  # target_include_directories(magnumafcpp PUBLIC "${PROJECT_BINARY_DIR}" )

  # Must use GNUInstallDirs to install libraries into correct locations on all
  # platforms.
  include(GNUInstallDirs)

  # 'make install' to the correct locations (provided by GNUInstallDirs).
  install(
    TARGETS magnumafcpp
    EXPORT magnumafcppConfig
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # This is for Windows
  )

  install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.hpp")
  install(DIRECTORY ${PROJECT_BINARY_DIR}/include/
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}) # this includes cmake
                                                   # generated headers e.g. the
                                                   # version header.

  # This makes the project importable from the install directory Put config file
  # in per-project dir (name MUST match), can also just go into 'cmake'.
  # install(EXPORT magnumafcppConfig DESTINATION share/magnumaf/cmake)
  install(EXPORT magnumafcppConfig DESTINATION cmake/)

  # This makes the project importable from the build directory
  export(TARGETS magnumafcpp FILE magnumafcppConfig.cmake)

  # Building Installer from:
  # https://cmake.org/cmake/help/latest/guide/tutorial/index.html#building-an-installer-step-7
  include(InstallRequiredSystemLibraries)
  set(CPACK_PACKAGE_VENDOR "Paul Heistracher")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
      "magnum.af: A finite differences GPU-accelerated micromagnetic and atomistic simulation software"
  )
  set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
  set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
  set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
  set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
  if(PACK_AS_DEB) # pack as .deb file, instead of package script and .tar files.
    set(CPACK_GENERATOR "DEB")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Paul Heistracher")
    set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON) # resolve dependencies automatically
    # WORKS#SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libvtk7.1 (>= 7.1.1), libvtk7.1 (<
    # 7.2)")
  endif()
  include(CPack)
endif()
