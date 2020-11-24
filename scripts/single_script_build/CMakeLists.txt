find_package(ArrayFire REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

file(
  GLOB EXECUTABLE
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.cpp)

# Build the program, linking specifically with designated backends ArrayFire CPU
# backend
if(${ArrayFire_CPU_FOUND})
  message(STATUS "ArrayFire CPU backend found. Enabling CPU benchmark")
  add_executable(magnum.af-cpu ${EXECUTABLE})
  target_link_libraries(magnum.af-cpu magnumafcpp ${ArrayFire_CPU_LIBRARIES}
                        ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} stdc++fs)
  # when linking agsinst lib use (linking here is slow, however):
  # ADD_EXECUTABLE(magnum.af-cpu ${EXECUTABLE})
  # TARGET_LINK_LIBRARIES(magnum.af-cpu magnumafcpp ${ArrayFire_CPU_LIBRARIES}
  # ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} stdc++fs)
endif()

# ArrayFire OpenCL backend
find_package(OpenCL)
if(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
  # We need to find OpenCL as transitive linking is disabled on some OSes
  message(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL benchmark")
  add_executable(magnum.af-opencl ${EXECUTABLE})
  target_link_libraries(
    magnum.af-opencl magnumafcpp ${ArrayFire_OpenCL_LIBRARIES} ${OpenCL_LIBRARIES}
    ${VTK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} stdc++fs)
endif()

# ArrayFire CUDA backend
find_package(CUDA)
if(${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
  # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
  find_package(CUDA REQUIRED)
  find_package(NVVM REQUIRED)
  message(STATUS ${CUDA_TOOLKIT_ROOT_DIR})
  message(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")
  add_executable(magnum.af-cuda ${EXECUTABLE})
  target_link_libraries(
    magnum.af-cuda
    magnumafcpp
    ${ArrayFire_CUDA_LIBRARIES}
    ${VTK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CUDA_LIBRARIES}
    ${NVVM_LIB}
    stdc++fs)
endif()
