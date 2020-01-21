CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

project(magnum.af VERSION 0.1.0) # sets magnum.af_VERSION_XXX and PROJECT_VERSION_XXX
MESSAGE( STATUS "Setting up project ${PROJECT_NAME} version ${PROJECT_VERSION}." )
#MESSAGE( STATUS "-------->>>> ${PROJECT_VERSION_MAJOR}" )
#MESSAGE( STATUS "-------->>>> ${PROJECT_VERSION_MINOR}" )
#MESSAGE( STATUS "-------->>>> ${PROJECT_VERSION_PATCH}" )
#MESSAGE( STATUS "-------->>>> ${PROJECT_VERSION_TWEAK}" )

option(EXAMPLES "Build examples." ON)
option(TESTS "Build unit and integration tests." ON)
option(PYTHON "Build python bindings" ON)
option(DEBUG "Build with debugging flag '-g'." OFF)
option(SETINSTALL "Building shared library, configuring install directories." ON)
option(PACK_AS_DEB "Pack as .deb file instead of package script and .tar files." ON)
option(CUDA "Build CUDA backend." ON)
option(OPENCL "Build OPENCL backend." ON)

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-copy") # silence implicitly-declared ‘af::dim4& af::dim4::operator=(const af::dim4&)’ is deprecated
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

## build library
# Generate Version Config Header
configure_file(${PROJECT_SOURCE_DIR}/include/magnumafcpp/magnumafConfig.hpp.in ${PROJECT_BINARY_DIR}/include/magnumafcpp/magnumafConfig.hpp)

# Generate Git Config header
set(PRE_CONFIGURE_FILE  "${PROJECT_SOURCE_DIR}/include/magnumafcpp/magnumafConfig_git.hpp.in" )
set(POST_CONFIGURE_FILE "${PROJECT_BINARY_DIR}/include/magnumafcpp/magnumafConfig_git.hpp" )
include(${PROJECT_SOURCE_DIR}/CMakeModules/git_watcher.cmake)

include_directories( ${PROJECT_SOURCE_DIR}/include/magnumafcpp
                     ${PROJECT_SOURCE_DIR}/include/magnumafcpp/llg_terms
                     ${PROJECT_SOURCE_DIR}/include/magnumafcpp/solvers
                     ${PROJECT_SOURCE_DIR}/include/magnumafcpp/util
                     ${PROJECT_SOURCE_DIR}/include/magnumafcpp/integrators
                     ${PROJECT_BINARY_DIR}/include/magnumafcpp/ # for configured headers
)

add_subdirectory(src)

# build examples
if(EXAMPLES)
    message(STATUS "Building example simulation scripts (option 'EXAMPLES').")
    add_subdirectory(examples)
endif()

# build tests
if(TESTS)
    message(STATUS "Enable unit and integration testing (option 'TESTING').")
    enable_testing() # enables CTest with 'make test' and 'ctest' commands
    add_subdirectory(test/unit) # adding unit tests
    add_subdirectory(test/integration) # adding integration tests
endif()

# build single scripts for simulations
add_subdirectory(scripts/single_script_build)

# build python bindings
if(PYTHON)
    message(STATUS "Building python bindings (option 'PYTHON_BINDINGS').")
    add_subdirectory(python)
endif()

# uninstall target
# alternatively, use shell for
# removing installed files:      $ cat install_manifest.txt | sudo xargs rm
# removing remaining empty dirs: $ cat install_manifest.txt | xargs -L1 dirname | sudo xargs rmdir -p --ignore-fail-on-non-empty 2>/dev/null
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake  COMMAND ${CMAKE_COMMAND} -E echo "Note: to remove remaining empty folders run: \\'cat install_manifest.txt \\| xargs -L1 dirname \\| sudo xargs rmdir -p --ignore-fail-on-non-empty 2\\>\\/dev\\/null\\'")
endif()
