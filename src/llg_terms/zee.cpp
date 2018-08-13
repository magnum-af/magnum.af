#include "zee.hpp"

Zee::Zee(af::array zee_in) : zee_field(zee_in)
{
}

Zee::Zee(long int aptr)
{
    void **a= (void**)aptr;
    zee_field = *( new af::array( *a ));
}

Zee::Zee(af::array (*callback_func_in)(State state)): callback_func(callback_func_in)
{
    callback=true;
}

af::array Zee::h(const State& state){
    //timer = timer::start();
    //if(param.afsync) sync();
    //time += timer::stop(timer);
    if(callback){
        return callback_func(state);
    }
    else return zee_field;
}

//Zeeman energy term
double Zee::E(const State& state){
    return - state.param.mu0 * state.param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy *state.mesh.dz; 
}
