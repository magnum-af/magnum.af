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

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Link runTests with what we want to test and the GTest and pthread library
file( GLOB APP_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp )
foreach( testsourcefile ${APP_SOURCES} )
    # String replace, to cut off .cpp.
    string( REPLACE ".cpp" "" testname ${testsourcefile} )
    # Adding executables and linking
    add_executable( ${testname} ${testsourcefile} )
    target_link_libraries( ${testname} ${GTEST_LIBRARIES} ${ArrayFire_CPU_LIBRARIES} ${VTK_LIBRARIES} pthread )
endforeach( testsourcefile ${APP_SOURCES} )
