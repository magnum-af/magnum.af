CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(magnum.af)

FIND_PACKAGE(ArrayFire REQUIRED)
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FILE(GLOB EXAMPLES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp )

# Build the program, linking specifically with designated backends

# ArrayFire CPU backend
foreach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
    # String replace, to get filename without extension
    string( REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE} )
    get_filename_component(EXAMPLE ${EXAMPLE} NAME )
    set(EXAMPLE ${EXAMPLE}_cpu)
    add_executable( ${EXAMPLE} ${EXAMPLE_SOURCE_FILE} )
    TARGET_LINK_LIBRARIES( ${EXAMPLE} magnum.af_CORE_SHARED ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endforeach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )

# ArrayFire OpenCL backend
FIND_PACKAGE(OpenCL)
IF(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
    # We need to find OpenCL as transitive linking is disabled on some OSes
    MESSAGE(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL benchmark")
    foreach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
        # String replace, to get filename without extension
        string( REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE} )
        get_filename_component(EXAMPLE ${EXAMPLE} NAME )
        set(EXAMPLE ${EXAMPLE}_opencl)
        add_executable( ${EXAMPLE} ${EXAMPLE_SOURCE_FILE} )
        TARGET_LINK_LIBRARIES(${EXAMPLE} magnum.af_CORE_SHARED ${ArrayFire_OpenCL_LIBRARIES} ${OpenCL_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endforeach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
ENDIF()

# ArrayFire CUDA backend
FIND_PACKAGE(CUDA)
IF(${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
    # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
    FIND_PACKAGE(CUDA REQUIRED)
    FIND_PACKAGE(NVVM REQUIRED)
    MESSAGE(STATUS ${CUDA_TOOLKIT_ROOT_DIR})
    MESSAGE(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")
    foreach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
        # String replace, to get filename without extension
        string( REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE} )
        get_filename_component(EXAMPLE ${EXAMPLE} NAME )
        set(EXAMPLE ${EXAMPLE}_TODO)
        add_executable( ${EXAMPLE} ${EXAMPLE_SOURCE_FILE} )
        TARGET_LINK_LIBRARIES(${EXAMPLE} magnum.af_CORE_SHARED ${ArrayFire_CUDA_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES} ${NVVM_LIB})
    endforeach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
ENDIF()
