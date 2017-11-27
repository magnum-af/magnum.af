#include "zee.hpp"
#include "func.hpp"

Zee::Zee(af::array zee_in, Mesh mesh_in, Param param_in) :
     zee_field(zee_in), mesh(mesh_in), param(param_in){
}

Zee::Zee(long int aptr, Mesh mesh_in, Param param_in): mesh(mesh_in), param(param_in){
    void **a= (void**)aptr;
    zee_field = *( new af::array( *a ));
}

array Zee::h(const State& state){
  //timer = timer::start();
  //if(param.afsync) sync();
  //time += timer::stop(timer);
  return zee_field;
}

//Zeeman energy term
double Zee::E(const State& state){
    return - param.mu0 * param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz; 
}
