#include "zee.hpp"
#include "func.hpp"

Zee::Zee(af::array zee_in) : zee_field(zee_in)
{
}

Zee::Zee(double rate_in, double hzee_max_in) : rate(rate_in), hzee_max(hzee_max_in)
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

array Zee::h(const State& state){
    //timer = timer::start();
    //if(param.afsync) sync();
    //time += timer::stop(timer);
    if(callback){
        return callback_func(state);
    }
    else if (zee_field.isempty()){
        //HACK
        double field_Tesla = 0;
        if(state.t < hzee_max/rate) field_Tesla = rate *state.t; 
        else if (state.t < 3*hzee_max/rate) field_Tesla = -rate *state.t + 2*hzee_max; 
	else if(state.t < 4*hzee_max/rate) field_Tesla = rate*state.t - 4*hzee_max; 
        else {field_Tesla = 0; std::cout << "WARNING ZEE time out of range" << std::endl;}
        array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
        zee(span,span,span,0)=constant(field_Tesla/state.param.mu0 ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        return  zee;
    }
    else return zee_field;
}

//Zeeman energy term
double Zee::E(const State& state){
    return - state.param.mu0 * state.param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy *state.mesh.dz; 
}
