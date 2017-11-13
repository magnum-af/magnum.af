#include "teststate.hpp"
testState::testState (Mesh mesh_in, Param param_in):
              mesh(mesh_in),param(param_in)
{
}

void testState::printn0(){
  std::cout<<"mesh.n0="<<mesh.n0  <<std::endl;
}
