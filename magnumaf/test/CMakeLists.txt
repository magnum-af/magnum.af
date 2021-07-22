find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Locate Arrayfire
find_package(ArrayFire REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Link runTests with what we want to test and the GTest and pthread library
file(
  GLOB APP_SOURCES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "*.cpp" "*/*.cpp")
message(STATUS "Source files found are: ${APP_SOURCES}")

foreach(testsourcefile ${APP_SOURCES})
  # String replace, to cut off .cpp.
  string(REPLACE ".cpp" "" testname ${testsourcefile})
  get_filename_component(testname ${testname} NAME)
  set(testname test_${testname}_cpu)
  # Adding executables and linking
  add_executable(${testname} ${testsourcefile})
  target_link_libraries(${testname} magnumafcpp ${GTEST_LIBRARIES}
                        ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES} pthread)
  add_test(
    NAME ${testname}
    COMMAND ${testname})
endforeach(testsourcefile ${APP_SOURCES})

find_package(OpenCL)
if(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
  message(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL tests")
  foreach(testsourcefile ${APP_SOURCES})
    # String replace, to cut off .cpp.
    string(REPLACE ".cpp" "" testname ${testsourcefile})
    get_filename_component(testname ${testname} NAME)
    set(testname test_${testname}_opencl)
    # Adding executables and linking
    add_executable(${testname} ${testsourcefile})
    target_link_libraries(
      ${testname} magnumafcpp ${GTEST_LIBRARIES} ${ArrayFire_OpenCL_LIBRARIES}
      ${OpenCL_LIBRARIES} ${VTK_LIBRARIES} pthread)
    add_test(
      NAME ${testname}
      COMMAND ${testname})
    # add_test(test_opencl_${testname} test_opencl_${testname})
  endforeach(testsourcefile ${APP_SOURCES})
endif()

find_package(CUDA)
if(${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
  # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
  find_package(CUDA REQUIRED)
  find_package(NVVM REQUIRED)
  message(STATUS "ArrayFire CUDA found. Enabling CUDA tests")
  foreach(testsourcefile ${APP_SOURCES})
    # String replace, to cut off .cpp.
    string(REPLACE ".cpp" "" testname ${testsourcefile})
    get_filename_component(testname ${testname} NAME)
    set(testname test_${testname}_cuda)
    # Adding executables and linking
    add_executable(${testname} ${testsourcefile})
    target_link_libraries(
      ${testname}
      magnumafcpp
      ${GTEST_LIBRARIES}
      ${ArrayFire_CUDA_LIBRARIES}
      ${CMAKE_THREAD_LIBS_INIT}
      ${CUDA_LIBRARIES}
      ${NVVM_LIB}
      ${VTK_LIBRARIES}
      pthread)
    add_test(
      NAME ${testname}
      COMMAND ${testname})
  endforeach(testsourcefile ${APP_SOURCES})
endif()
