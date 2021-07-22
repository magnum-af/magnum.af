# Writing compile_commands.json for youcompleteme (YCM) vim-plugin Note: fully
# worked only after removing all global and local .ycm_extra_conf.py files with
# .json comile file only in project dir Found this hint in
# https://www.reddit.com/r/emacs/comments/7zxl37/trouble_setting_up_ycmd/
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
  execute_process(
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
      ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)
endif()

# Set c++ compiler version
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(COMPILER_SUPPORTS_CXX17)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
else()
  message(
    STATUS
      "Error: The compiler ${CMAKE_CXX_COMPILER} has no c++17 support. Please use a different C++ compiler or update the existing one."
  )
endif()

# avoid redundant linking when modifying shared libraries
# (especially when an executables source file changes)
set(CMAKE_LINK_DEPENDS_NO_SHARED TRUE)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -M") # test linking
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -H") # test linking tree
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpessimizing-move")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wredundant-move")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-conversion")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow") # shadowing caused a nasty bug
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")

if(DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()


# indicate we wish for the binaries to be written to the bin directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
# we have custom CMake configuration files, tell CMake to use them:
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

# build library Generate Version Config Header
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/magnumafcpp/magnumafConfig.hpp.in
               ${CMAKE_BINARY_DIR}/include/magnumafcpp/magnumafConfig.hpp)

# Generate Git Config header
set(PRE_CONFIGURE_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/include/magnumafcpp/magnumafConfig_git.hpp.in")
set(POST_CONFIGURE_FILE
    "${CMAKE_BINARY_DIR}/include/magnumafcpp/magnumafConfig_git.hpp")
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/git_watcher.cmake)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/magnumafcpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/magnumafcpp/field_terms
  ${CMAKE_CURRENT_SOURCE_DIR}/include/magnumafcpp/solvers
  ${CMAKE_CURRENT_SOURCE_DIR}/include/magnumafcpp/util
  ${CMAKE_CURRENT_SOURCE_DIR}/include/magnumafcpp/integrators
  ${CMAKE_BINARY_DIR}/include/magnumafcpp/ # for configured headers
)

add_subdirectory(src)

# build examples
if(EXAMPLES)
  message(STATUS "Building example simulation scripts (option 'EXAMPLES').")
  add_subdirectory(examples)
endif()

# build tests
if(TESTS)
  add_subdirectory(test) # adding unit tests
endif()

# uninstall target alternatively, use shell for removing installed files:      $
# cat install_manifest.txt | sudo xargs rm removing remaining empty dirs: $ cat
# install_manifest.txt | xargs -L1 dirname | sudo xargs rmdir -p
# --ignore-fail-on-non-empty 2>/dev/null
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)

  add_custom_target(
    uninstall
    COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    COMMAND
      ${CMAKE_COMMAND} -E echo
      "Note: to remove remaining empty folders run: \\'cat install_manifest.txt \\| xargs -L1 dirname \\| sudo xargs rmdir -p --ignore-fail-on-non-empty 2\\>\\/dev\\/null\\'"
  )
endif()
