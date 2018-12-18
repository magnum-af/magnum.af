#include "zee.hpp"

Zee::Zee(af::array zee_in) : zee_field(zee_in)
{
}

Zee::Zee(long int aptr)
{
    void **a= (void**)aptr;
    zee_field = *( new af::array( *a ));
}

long int Zee::get_m_addr(){
    af::array *a = new af::array(zee_field);
    return (long int) a->get();
}

Zee::Zee(af::array (*callback_func_in)(State state)): callback_func(callback_func_in)
{
    callback=true;
}

Zee::Zee(std::function<af::array(State)> lamda_callback): lamda_callback(lamda_callback){
    
    is_lamda=true;

}

///< Sets internal af::array to a global field (x, y, z) for all spacial dimensions given by state.mesh.dims().
void Zee::set_xyz(const State& state, const double x, const double y, const double z){
    af::dim4 dim = af::dim4(state.mesh.n0, state.mesh.n1, state.mesh.n2, 1);
    af::array temp = af::array(state.mesh.dims, f64);
    temp (af::span, af::span, af::span, 0) = af::constant(x, dim, f64);
    temp (af::span, af::span, af::span, 1) = af::constant(y, dim, f64);
    temp (af::span, af::span, af::span, 2) = af::constant(z, dim, f64);
    zee_field = temp;
}

af::array Zee::h(const State& state){
    //timer = timer::start();
    //if(param.afsync) sync();
    //time += timer::stop(timer);
    if (is_lamda) {
        return lamda_callback(state); 
    }
    else if(callback){
        return callback_func(state);
    }
    else return zee_field;
}

//Zeeman energy term
double Zee::E(const State& state){
    return - state.param.mu0 * state.param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy *state.mesh.dz;
}

double Zee::E(const State& state, const af::array& h){
    return - state.param.mu0 * state.param.ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy *state.mesh.dz;
}
