#include "state.hpp"
State::State (Mesh mesh_in, Param param_in, af::array m_in):
              mesh(mesh_in),param(param_in), m(m_in)
{
}

//State::State& operator+(af::array& rhs){
//  this->m=this->m+rhs;
//  return *this;
//}
