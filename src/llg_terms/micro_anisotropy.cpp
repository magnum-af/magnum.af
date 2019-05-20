//Ref Diss Abert p 15 sec 2.3 eq 2.3.28
//Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m) 
//With higher order (not implemented): Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m) ( + 4 Ku2 /(mu0 Ms) eu ( eu . m)^3
#include "micro_anisotropy.hpp"

//Energy calculation
//Edemag=-mu0/2 integral(M . Hdemag) dx
double UniaxialAnisotropyField::E(const State& state){
    return -constants::mu0/2. * state.material.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}

double UniaxialAnisotropyField::E(const State& state, const af::array& h){
    return -constants::mu0/2. * state.material.ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}


UniaxialAnisotropyField::UniaxialAnisotropyField (double Ku1, std::array<double, 3> Ku1_axis) : Ku1(Ku1), Ku1_axis(get_normalized_vector(Ku1_axis)){
}


UniaxialAnisotropyField::UniaxialAnisotropyField (af::array Ku1_field, std::array<double, 3> Ku1_axis) : Ku1_field(Ku1_field), Ku1_axis(get_normalized_vector(Ku1_axis)){
}


// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField (double Ku1, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2) : Ku1(Ku1), Ku1_axis(get_normalized_vector(std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2})){
}


// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField (long int Ku1_field_ptr, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2) : Ku1_field(*( new af::array( *((void**) Ku1_field_ptr)))), Ku1_axis(get_normalized_vector(std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2})){
}

af::array UniaxialAnisotropyField::h(const State& state){
    return calc_heff(state);
}

long int UniaxialAnisotropyField::h_ptr(const State& state){
    return (long int) (new af::array(calc_heff(state)))->get();
}

af::array UniaxialAnisotropyField::calc_heff(const State& state){
    af::timer timer_anisotropy = af::timer::start();

    af::array eu = af::array(state.mesh.dims, f64);//Normal vector
    eu(af::span,af::span,af::span,0) = Ku1_axis[0];
    eu(af::span,af::span,af::span,1) = Ku1_axis[1];
    eu(af::span,af::span,af::span,2) = Ku1_axis[2];

    af::array anisotropy = eu*state.m;
    anisotropy=af::sum(anisotropy,3);
    anisotropy=af::tile(anisotropy,1,1,1,3);

    if(state.afsync) af::sync();
    computation_time_heff += af::timer::stop(timer_anisotropy);
    if (state.Ms_field.isempty() && Ku1_field.isempty()){
        return  2.* Ku1/(constants::mu0 * state.material.ms) * (eu* anisotropy);
    }
    else if ( ! state.Ms_field.isempty() && Ku1_field.isempty()){
        af::array result =  2.* Ku1/(constants::mu0 * state.Ms_field) * (eu* anisotropy);
        af::replace(result, !af::iszero(state.Ms_field), 0); // Replacing all resulting NaN with 0
        return result;
    }
    else if ( state.Ms_field.isempty() && ! Ku1_field.isempty()){
        return  2.* Ku1_field/(constants::mu0 * state.material.ms) * (eu* anisotropy);
    }
    else {
        af::array result =  2.* Ku1_field/(constants::mu0 * state.Ms_field) * (eu* anisotropy);
        af::replace(result, !af::iszero(state.Ms_field), 0); // Replacing all resulting NaN with 0
        return result;
    }
}


std::array<double, 3> UniaxialAnisotropyField::get_normalized_vector(std::array<double, 3> vector){
    double norm = sqrt(pow(vector[0],2)+ pow(vector[1],2) + pow(vector[2], 2));
    return std::array<double, 3> {vector[0]/norm, vector[1]/norm, vector[2]/norm};
}


double UniaxialAnisotropyField::get_ku1_axis(int i){
    return Ku1_axis[i];
}


//void UniaxialAnisotropyField::set_Ku1_field(long int aptr){
//    void **a = (void **)aptr;
//    micro_Ku1_field = *( new af::array( *a ));
//}


long int UniaxialAnisotropyField::get_Ku1_field(){
    af::array *a = new af::array(Ku1_field);
    return (long int) a->get();
}
