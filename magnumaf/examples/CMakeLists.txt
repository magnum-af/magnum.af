find_package(ArrayFire REQUIRED)
find_package(VTK REQUIRED)
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

include(${VTK_USE_FILE})

file(
  GLOB_RECURSE EXAMPLES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.cpp)

# Linking with ArrayFire unified backend
foreach(EXAMPLE_SOURCE_FILE ${EXAMPLES})
  string(REPLACE ".cpp" "" EXAMPLE ${EXAMPLE_SOURCE_FILE})
  # stripping off possible paths to get a valid target name
  get_filename_component(EXAMPLE ${EXAMPLE} NAME)
  # MESSAGE("${EXAMPLE}")
  add_executable(${EXAMPLE} ${EXAMPLE_SOURCE_FILE})
  # Note: ArrayFire::af would only work when arrayfire unified is installed, but not e.g. cpu only:
  target_link_libraries(${EXAMPLE} ${ArrayFire_LIBRARIES} magnumafcpp ${VTK_LIBRARIES}
                        ${Boost_LIBRARIES})
endforeach(EXAMPLE_SOURCE_FILE ${EXAMPLES})

# Adding magnumafcpp_info binary to installation
install(TARGETS magnumafcpp_info RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
