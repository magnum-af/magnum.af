CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(magnum.af)# needed for 'inheritance'

# If the pyx file is a C++ file, we should specify that here.
set_source_files_properties(
  ${magnum.af_SOURCE_DIR}/magnumaf.pyx
  # MESSAGE(STATUS "----> PTH: PATH : ${magnum.af_SOURCE_DIR}")
  #Additional files would come here as follows:
  #${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/FurtherFiles.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )

cython_add_module_w_AF_VTK( magnumaf magnumaf.pyx $<TARGET_OBJECTS:magnum.af_CORE>)
# Further modules would come with cython_add_module_w_AF_VTK(target target.pyx [sources])
