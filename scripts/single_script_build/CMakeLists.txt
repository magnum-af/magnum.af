CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(magnum.af)

FIND_PACKAGE(ArrayFire REQUIRED)
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FILE(GLOB EXECUTABLE RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp )

# Build the program, linking specifically with designated backends
# ArrayFire CPU backend
if(${ArrayFire_CPU_FOUND})
    MESSAGE(STATUS "ArrayFire CPU backend found. Enabling CPU benchmark")
    ADD_EXECUTABLE(magnum.af-cpu ${EXECUTABLE})
    TARGET_LINK_LIBRARIES(magnum.af-cpu magnum.af_CORE_SHARED ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

# ArrayFire OpenCL backend
FIND_PACKAGE(OpenCL)
IF(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
    # We need to find OpenCL as transitive linking is disabled on some OSes
    MESSAGE(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL benchmark")
    ADD_EXECUTABLE(magnum.af-opencl  ${EXECUTABLE})
    TARGET_LINK_LIBRARIES(magnum.af-opencl magnum.af_CORE_SHARED ${ArrayFire_OpenCL_LIBRARIES} ${OpenCL_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

# ArrayFire CUDA backend
FIND_PACKAGE(CUDA)
IF(${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
    # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
    FIND_PACKAGE(CUDA REQUIRED)
    FIND_PACKAGE(NVVM REQUIRED)
    MESSAGE(STATUS ${CUDA_TOOLKIT_ROOT_DIR})
    MESSAGE(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")
    ADD_EXECUTABLE(magnum.af-cuda  ${EXECUTABLE})
    TARGET_LINK_LIBRARIES(magnum.af-cuda magnum.af_CORE_SHARED ${ArrayFire_CUDA_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES} ${NVVM_LIB})
ENDIF()
