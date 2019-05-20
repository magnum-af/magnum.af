#include "zee.hpp"

ExternalField::ExternalField(af::array zee_in) : zee_field(zee_in)
{
}

ExternalField::ExternalField(long int aptr) : zee_field(*( new af::array( *((void**) aptr))))
{
}

long int ExternalField::get_m_addr(){
    af::array *a = new af::array(zee_field);
    return (long int) a->get();
}

ExternalField::ExternalField(af::array (*callback_func_in)(State state)): callback_func(callback_func_in), callback(true)
{
}

ExternalField::ExternalField(std::function<af::array(State)> lamda_callback): lamda_callback(lamda_callback), is_lamda(true)
{
}

///< Sets internal af::array to a global field (x, y, z) for all spacial dimensions given by state.mesh.dims().
void ExternalField::set_homogenuous_field(const double x, const double y, const double z){
    af::dim4 dim = af::dim4(zee_field.dims(0), zee_field.dims(1), zee_field.dims(2), 1);
    zee_field(af::span, af::span, af::span, 0) = af::constant(x, dim, f64);
    zee_field(af::span, af::span, af::span, 1) = af::constant(y, dim, f64);
    zee_field(af::span, af::span, af::span, 2) = af::constant(z, dim, f64);
}

af::array ExternalField::h(const State& state){
    //af::timer timer = af::timer::start();
    //if(state.afsync) sync();
    //af_time += af::timer::stop(timer);
    if (is_lamda) {
        return lamda_callback(state); 
    }
    else if(callback){
        return callback_func(state);
    }
    else return zee_field;
}

//Zeeman energy term
double ExternalField::E(const State& state){
    return - constants::mu0 * state.material.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy *state.mesh.dz;
}

double ExternalField::E(const State& state, const af::array& h){
    return - constants::mu0 * state.material.ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy *state.mesh.dz;
}
