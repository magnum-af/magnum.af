#include "state.hpp"
State::State (Mesh mesh_in, Param param_in, af::array m_in):
              mesh(mesh_in),param(param_in), m(m_in)
{
}

State::State (Mesh mesh_in, Param param_in, long int aptr):
              mesh(mesh_in),param(param_in)
{
  void **a = (void **)aptr;
  m = *( new af::array( *a ));
  m.lock();
}

void State::_vti_writer_micro(std::string outputname){
    vti_writer_micro(m, mesh, outputname); 
}
void State::_vti_writer_atom (std::string outputname){
    vti_writer_atom(m, mesh, outputname); 
}
void State::_vti_reader(std::string inputname){
    vti_reader(m, mesh, inputname);
}


void State::_vtr_writer(std::string outputname){
    vtr_writer(m, mesh, outputname); 
}
void State::_vtr_reader(std::string inputname){
    vtr_reader(m, mesh, inputname); 
}
