find_package(ArrayFire REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

file(
  GLOB EXAMPLES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.cpp)

# Build the program, linking specifically with designated backends

# ArrayFire CPU backend
foreach(EXAMPLE_SOURCE_FILE ${EXAMPLES})
  # String replace, to get filename without extension
  string(REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE})
  get_filename_component(EXAMPLE ${EXAMPLE} NAME)
  set(EXAMPLE ${EXAMPLE}_cpu)
  add_executable(${EXAMPLE} ${EXAMPLE_SOURCE_FILE}
                            $<TARGET_OBJECTS:magnum.af_CORE_OBJECT>)
  target_link_libraries(
    ${EXAMPLE} ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}) # add magnum.af_CORE_SHARED here  when moving
                               # from OBJECT to SHARED lib
endforeach(EXAMPLE_SOURCE_FILE ${EXAMPLES})

# ArrayFire OpenCL backend
find_package(OpenCL)
if(${OPENCL}
   AND ${ArrayFire_OpenCL_FOUND}
   AND ${OpenCL_FOUND})
  # We need to find OpenCL as transitive linking is disabled on some OSes
  message(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL benchmark")
  foreach(EXAMPLE_SOURCE_FILE ${EXAMPLES})
    # String replace, to get filename without extension
    string(REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE})
    get_filename_component(EXAMPLE ${EXAMPLE} NAME)
    set(EXAMPLE ${EXAMPLE}_opencl)
    add_executable(${EXAMPLE} ${EXAMPLE_SOURCE_FILE}
                              $<TARGET_OBJECTS:magnum.af_CORE_OBJECT>)
    target_link_libraries(
      ${EXAMPLE} ${ArrayFire_OpenCL_LIBRARIES} ${OpenCL_LIBRARIES}
      ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  endforeach(EXAMPLE_SOURCE_FILE ${EXAMPLES})
endif()

# ArrayFire CUDA backend
find_package(CUDA)
if(${CUDA}
   AND ${ArrayFire_CUDA_FOUND}
   AND ${CUDA_FOUND})
  # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
  find_package(CUDA REQUIRED)
  find_package(NVVM REQUIRED)
  message(STATUS ${CUDA_TOOLKIT_ROOT_DIR})
  message(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")
  foreach(EXAMPLE_SOURCE_FILE ${EXAMPLES})
    # String replace, to get filename without extension
    string(REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE})
    get_filename_component(EXAMPLE ${EXAMPLE} NAME)
    set(EXAMPLE ${EXAMPLE}_cuda)
    add_executable(${EXAMPLE} ${EXAMPLE_SOURCE_FILE}
                              $<TARGET_OBJECTS:magnum.af_CORE_OBJECT>)
    target_link_libraries(
      ${EXAMPLE} ${ArrayFire_CUDA_LIBRARIES} ${VTK_LIBRARIES}
      ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES} ${NVVM_LIB})
  endforeach(EXAMPLE_SOURCE_FILE ${EXAMPLES})
endif()

# Adding magnumafcpp_info binary to installation
install(TARGETS magnumafcpp_info_cpu
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(${OPENCL}
   AND ${ArrayFire_OpenCL_FOUND}
   AND ${OpenCL_FOUND})
  install(TARGETS magnumafcpp_info_opencl
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
if(${CUDA}
   AND ${ArrayFire_CUDA_FOUND}
   AND ${CUDA_FOUND})
  install(TARGETS magnumafcpp_info_cuda
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
