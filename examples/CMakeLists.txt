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
    add_executable( ${EXAMPLE} ${EXAMPLE_SOURCE_FILE} $<TARGET_OBJECTS:magnum.af_CORE_OBJECT>)
    TARGET_LINK_LIBRARIES( ${EXAMPLE} ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}) # add magnum.af_CORE_SHARED here  when moving from OBJECT to SHARED lib
endforeach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )

# ArrayFire OpenCL backend
FIND_PACKAGE(OpenCL)
IF(${OPENCL} AND ${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
    # We need to find OpenCL as transitive linking is disabled on some OSes
    MESSAGE(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL benchmark")
    foreach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
        # String replace, to get filename without extension
        string( REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE} )
        get_filename_component(EXAMPLE ${EXAMPLE} NAME )
        set(EXAMPLE ${EXAMPLE}_opencl)
        add_executable( ${EXAMPLE} ${EXAMPLE_SOURCE_FILE} $<TARGET_OBJECTS:magnum.af_CORE_OBJECT>)
        TARGET_LINK_LIBRARIES(${EXAMPLE} ${ArrayFire_OpenCL_LIBRARIES} ${OpenCL_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endforeach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
ENDIF()

# ArrayFire CUDA backend
FIND_PACKAGE(CUDA)
IF(${CUDA} AND ${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
    # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
    FIND_PACKAGE(CUDA REQUIRED)
    FIND_PACKAGE(NVVM REQUIRED)
    MESSAGE(STATUS ${CUDA_TOOLKIT_ROOT_DIR})
    MESSAGE(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")
    foreach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
        # String replace, to get filename without extension
        string( REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE} )
        get_filename_component(EXAMPLE ${EXAMPLE} NAME )
        set(EXAMPLE ${EXAMPLE}_cuda)
        add_executable( ${EXAMPLE} ${EXAMPLE_SOURCE_FILE} $<TARGET_OBJECTS:magnum.af_CORE_OBJECT>)
        TARGET_LINK_LIBRARIES(${EXAMPLE} ${ArrayFire_CUDA_LIBRARIES} ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES} ${NVVM_LIB})
    endforeach( EXAMPLE_SOURCE_FILE ${EXAMPLES} )
ENDIF()

# Adding magnumafcpp_info binary to installation
install(TARGETS magnumafcpp_info_cpu RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
if( ${OPENCL} )
install(TARGETS magnumafcpp_info_opencl RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif()
if( ${CUDA} )
install(TARGETS magnumafcpp_info_cuda RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif()
