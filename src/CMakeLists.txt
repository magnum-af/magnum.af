CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(pth-mag)

#useage:
#on GTO use:
#cmake -DVTK_DIR:PATH=/home/paul/Programs/VKT-build ..

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
#set(VTK_DIR "/usr/include/vtk-6.2/")#TODO Unluckily this has to be set for each system manually, according to VTK
#set(VTK_DIR "/usr/include/vtk-5.8/")# for GTO #TODO Unluckily this has to be set for each system manually, according to VTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})


# Include the ArrayFire hreaders
INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS})

# Get a list of all of the source files in the directory:
FILE(GLOB SOURCE "*.cpp")

# There are several ways of compiling source code in CMake. In most cases you
# specify the source files to an ADD_EXCUTABLE call. Because we intend to
# link against multiple ArrayFire backends, this approach would result in
# every source file being re-compiled multiple times.
# Here we take a different approach to compilation. We build all of the source
# files into a CMake object library and then link this library against the
# different ArrayFire backends.
ADD_LIBRARY(pth-mag_CORE OBJECT ${SOURCE})

# Build the program, linking specifically with designated backends
# ArrayFire CPU backend
if(${ArrayFire_CPU_FOUND})
    MESSAGE(STATUS "ArrayFire CPU backend found. Enabling CPU benchmark")
    ADD_EXECUTABLE(pth-mag-cpu $<TARGET_OBJECTS:pth-mag_CORE>)
    TARGET_LINK_LIBRARIES(pth-mag-cpu ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

# ArrayFire OpenCL backend
FIND_PACKAGE(OpenCL)
MESSAGE(STATUS "----> PTH: OpenCL_FOUND: ${OpenCL_FOUND}")
MESSAGE(STATUS "----> PTH:             : ${OpenCL_INCLUDE_DIRS}")
MESSAGE(STATUS "----> PTH:             : ${OpenCL_LIBRARIES}")
MESSAGE(STATUS "----> PTH:             : ${OpenCL_VERSION_STRING}")
MESSAGE(STATUS "----> PTH:             : ${OpenCL_VERSION_MAJOR}")
MESSAGE(STATUS "----> PTH:             : ${OpenCL_VERSION_MINOR}")
MESSAGE(STATUS "----> PTH:             : ${OpenCL_INCLUDE_DIR}")
MESSAGE(STATUS "----> PTH:             : ${OpenCL_LIBRARY}")
#IF(${ArrayFire_OpenCL_FOUND})
# TODO OpenCL not found by FIND_PACKAGE on this machine only but included by arrayfire
IF(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
    # We need to find OpenCL as transitive linking is disabled on some OSes
    MESSAGE(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL benchmark")
    ADD_EXECUTABLE(pth-mag-opencl $<TARGET_OBJECTS:pth-mag_CORE>)
    TARGET_LINK_LIBRARIES(pth-mag-opencl ${ArrayFire_OpenCL_LIBRARIES} 
        # TODO OpenCL not found by FIND_PACKAGE on this machine only but included by arrayfire
        ${OpenCL_LIBRARIES} #TODO include this on eg GTO
        ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

# ArrayFire CUDA backend
FIND_PACKAGE(CUDA)
IF(${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
    # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
    FIND_PACKAGE(CUDA REQUIRED)
    FIND_PACKAGE(NVVM REQUIRED)
    MESSAGE(STATUS ${CUDA_TOOLKIT_ROOT_DIR})
    MESSAGE(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")
    ADD_EXECUTABLE(pth-mag-cuda $<TARGET_OBJECTS:pth-mag_CORE>)
    TARGET_LINK_LIBRARIES(pth-mag-cuda ${ArrayFire_CUDA_LIBRARIES} ${VTK_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES} ${NVVM_LIB})
ENDIF()


# If the pyx file is a C++ file, we should specify that here.
set_source_files_properties(
  ${pth-mag_SOURCE_DIR}/pth_mag.pyx
  # MESSAGE(STATUS "----> PTH: PATH : ${pth-mag_SOURCE_DIR}")
  #Additional files would come here as follows:
  #${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/FurtherFiles.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )


cython_add_module_w_AF_VTK( pth_mag pth_mag.pyx param.cpp controller.cpp integrator.cpp mesh.cpp state.cpp micro_demag.cpp func.cpp llg.cpp micro_exch.cpp zee.cpp vtk_IO.cpp)
# Further modules would come with cython_add_module_w_AF_VTK(target target.pyx [sources])







    #if(VTK_LIBRARIES)
    #  target_link_libraries(pth-mag-opencl ${VTK_LIBRARIES})
    #else()
    #  target_link_libraries(pth-mag-opencl vtkHybrid vtkWidgets)
    #endif()
