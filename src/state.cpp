#include "state.hpp"
State::State (Mesh mesh_in, Param param_in, af::array m_in):
              mesh(mesh_in),param(param_in), m(m_in)
{
    m.lock();
}

State::State (Mesh mesh_in, Param param_in, long int aptr):
              mesh(mesh_in),param(param_in)
{
  void **a = (void **)aptr;
  m = *( new af::array( *a ));
  //af::array *A = new af::array(*a); //Alternative step-by-step
  //m=*A;			      //
  m.lock();
}

//Causes problem with wrapping
void State::write_vtk(std::string outputname){
   vtk_writer(m, mesh, outputname); 
}

void State::write_vti(std::string outputname){
   vti_writer_micro(m, mesh, outputname); 
}

void State::write_vtk_todel(){
   vtk_writer(m, mesh, "test"); 
}

//State::State& operator+(af::array& rhs){
//  this->m=this->m+rhs;
//  return *this;
//}
