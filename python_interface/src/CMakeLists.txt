# If the pyx file is a C++ file, we should specify that here.
set_source_files_properties(
  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/pth_mag.pyx
  #Additional files would come here
  #${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/FurtherFiles.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )


cython_add_module_w_AF_VTK( pth_mag pth_mag.pyx ../../src/mesh.cpp ../../src/state.cpp ../../src/micro_demag.cpp ../../src/func.cpp ../../src/llg.cpp ../../src/micro_exch.cpp ../../src/zee.cpp ../../src/vtk_writer.cpp)
# Further modules would come with cython_add_module_w_AF_VTK(target target.pyx [sources])

# Multi-file cython modules do not appear to be working at the moment.
#cython_add_module( rect rect.pyx Rectangle.cxx )
#INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS}) # When is this needed?
#If VTK_DIR not set, you may set: 
#set(VTK_DIR "/home/pheistra/Programs/VTK-build/")
#set( PYTHON_LIBRARIES ${PYTHON_LIBRARIES} ${ArrayFire_OpenCL_LIBRARIES} ${VTK_LIBRARIES})
#FIND_PACKAGE(OpenCL)
