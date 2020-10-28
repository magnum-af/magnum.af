# Compile with C++11 support.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(
    STATUS
      "Error: The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler or update the existing one."
  )
endif()

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Locate Arrayfire
find_package(ArrayFire REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Link runTests with what we want to test and the GTest and pthread library
file(
  GLOB APP_SOURCES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.cpp)
message(STATUS "Source files found are: ${APP_SOURCES}")

foreach(testsourcefile ${APP_SOURCES})
  # String replace, to cut off .cpp.
  string(REPLACE ".cpp" "" testname ${testsourcefile})
  get_filename_component(testname ${testname} NAME)
  set(testname test_cpu_${testname})
  # Adding executables and linking
  add_executable(${testname} ${testsourcefile}
                             $<TARGET_OBJECTS:magnumaf_object_lib>)
  target_link_libraries(${testname} ${GTEST_LIBRARIES}
                        ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES} pthread)
  add_test(
    NAME ${testname}
    COMMAND ${testname}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
endforeach(testsourcefile ${APP_SOURCES})

find_package(OpenCL)
if(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
  message(STATUS "ArrayFire OpenCL backend found. Enabling OpenCL tests")
  foreach(testsourcefile ${APP_SOURCES})
    # String replace, to cut off .cpp.
    string(REPLACE ".cpp" "" testname ${testsourcefile})
    get_filename_component(testname ${testname} NAME)
    set(testname test_opencl_${testname})
    # Adding executables and linking
    add_executable(${testname} ${testsourcefile}
                               $<TARGET_OBJECTS:magnumaf_object_lib>)
    target_link_libraries(
      ${testname} ${GTEST_LIBRARIES} ${ArrayFire_OpenCL_LIBRARIES}
      ${OpenCL_LIBRARIES} ${VTK_LIBRARIES} pthread)
    add_test(
      NAME ${testname}
      COMMAND ${testname}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
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
    set(testname test_cuda_${testname})
    # Adding executables and linking
    add_executable(${testname} ${testsourcefile}
                               $<TARGET_OBJECTS:magnumaf_object_lib>)
    target_link_libraries(
      ${testname}
      ${GTEST_LIBRARIES}
      ${ArrayFire_CUDA_LIBRARIES}
      ${CMAKE_THREAD_LIBS_INIT}
      ${CUDA_LIBRARIES}
      ${NVVM_LIB}
      ${VTK_LIBRARIES}
      pthread)
    add_test(
      NAME ${testname}
      COMMAND ${testname}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
  endforeach(testsourcefile ${APP_SOURCES})
endif()
