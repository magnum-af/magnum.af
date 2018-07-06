cmake_minimum_required(VERSION 2.6)

# Compile with C++11 support.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "Error: The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler or update the existing one.")
endif()
 
# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Locate Arrayfire
FIND_PACKAGE(ArrayFire REQUIRED)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Link runTests with what we want to test and the GTest and pthread library
add_executable(param tparam.cpp)
target_link_libraries(param ${GTEST_LIBRARIES} ${ArrayFire_CPU_LIBRARIES} pthread)

add_executable(mesh tmesh.cpp)
target_link_libraries(mesh ${GTEST_LIBRARIES} ${ArrayFire_CPU_LIBRARIES} pthread)
