#include "zee.hpp"
#include "func.hpp"

Zee::Zee(af::array zee_in) : zee_field(zee_in)
{
}

Zee::Zee(double rate_in) : rate(rate_in)
{
}

Zee::Zee(long int aptr)
{
    void **a= (void**)aptr;
    zee_field = *( new af::array( *a ));
}

array Zee::h(const State& state){
    //timer = timer::start();
    //if(param.afsync) sync();
    //time += timer::stop(timer);
    if (zee_field.isempty()){
        array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
        zee(span,span,span,0)=constant(rate/state.param.mu0 *state.t,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        return  zee;
    }
    else return zee_field;
}

//Zeeman energy term
double Zee::E(const State& state){
    return - state.param.mu0 * state.param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy *state.mesh.dz; 
}
