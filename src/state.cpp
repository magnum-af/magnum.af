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


void State::print_m(){
    af::print("m=", m);
}

//TODO
//long int State::get_m(){
//    return (long int)&m;
//    //return &m.arr;
//}
//TODO END

//State::State& operator+(af::array& rhs){
//  this->m=this->m+rhs;
//  return *this;
//}
