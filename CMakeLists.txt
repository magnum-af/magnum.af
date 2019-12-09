CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(magnum.af)

option(EXAMPLES "Build examples." ON)
option(TESTING "Build unit and integration tests." ON)
option(PYTHON_BINDINGS "Build python bindings" ON)
option(DEBUG "Build with debugging flag '-g'." OFF)

# Writing compile_commands.json for youcompleteme (YCM) vim-plugin
# Note: fully worked only after removing all global and local .ycm_extra_conf.py files with .json comile file only in project dir
# Found this hint in https://www.reddit.com/r/emacs/comments/7zxl37/trouble_setting_up_ycmd/
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

# Compile with C++11 support.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "Error: The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler or update the existing one.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmisleading-indentation")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wduplicated-cond")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wduplicated-branches")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wlogical-op")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnull-dereference")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdouble-promotion")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wformat=2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-align")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-aliasing")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
#TODO#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-conversion")
#TODO#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion")
#TODO mainly initializer list#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")

if(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# indicate we wish for the binaries to be written to the bin directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
# we have custom CMake configuration files, tell CMake to use them:
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# build library
add_subdirectory(src)

# build examples
if(EXAMPLES)
    message(STATUS "Building example simulation scripts (option 'EXAMPLES').")
    add_subdirectory(examples)
endif()

# build tests
if(TESTING)
    message(STATUS "Enable unit and integration testing (option 'TESTING').")
    enable_testing() # enables CTest with 'make test' and 'ctest' commands
    add_subdirectory(test/unit) # adding unit tests
    add_subdirectory(test/integration) # adding integration tests
endif()

# build single scripts for simulations
add_subdirectory(scripts/single_script_build)

# build python bindings
if(PYTHON_BINDINGS)
    message(STATUS "Building python bindings (option 'PYTHON_BINDINGS').")
    add_subdirectory(python)
endif()
