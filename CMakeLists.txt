CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(pth-mag)
#Note in new projects, set VTK_DIR by export VTK_DIR=/home/.../VTK-build:$VTK_DIR

# Compile with C++14 support.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
    message(STATUS "Error: The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler or update the existing one.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# indicate we wish for the binaries to be written to the bin directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
# we have custom CMake configuration files, tell CMake to use them:
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# build items in subdirectories
add_subdirectory(src)
