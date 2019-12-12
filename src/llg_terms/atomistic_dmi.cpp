#include "atomistic_dmi.hpp"
#include "../func.hpp"

namespace magnumafcpp{

double AtomisticDmiField::E(const State& state){
  return -constants::mu0/2. *state.Ms * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3));
}

double AtomisticDmiField::E(const State& state, const af::array& h){
  return -constants::mu0/2. *state.Ms * afvalue(sum(sum(sum(sum(h * state.m, 0), 1), 2), 3));
}


std::array<double, 3> get_normalized_vector(std::array<double, 3> vector){
    double norm = sqrt(pow(vector[0], 2)+ pow(vector[1], 2) + pow(vector[2], 2));
    return std::array<double, 3> {vector[0]/norm, vector[1]/norm, vector[2]/norm};
}


AtomisticDmiField::AtomisticDmiField (const double D_atom, std::array<double, 3> D_atom_axis) : D_atom(D_atom), D_atom_axis(get_normalized_vector(D_atom_axis))
{
}

af::array AtomisticDmiField::h(const State& state){
    af::timer timer_dmi = af::timer::start();

    af::array n= af::array(state.mesh.n0, state.mesh.n1, state.mesh.n2, 3, f64);
    n(af::span, af::span, af::span, 0)=D_atom_axis[0];
    n(af::span, af::span, af::span, 1)=D_atom_axis[1];
    n(af::span, af::span, af::span, 2)=D_atom_axis[2];
    //print("n", n);

    af::array filtr_fd1=af::constant(0.0, 3, 3, 3, 3, f64);
    //dmx/dx
    filtr_fd1(0, 1, 1, 0)= 1;
    filtr_fd1(2, 1, 1, 0)=-1;

    //dmy/dy
    filtr_fd1(1, 0, 1, 1)= 1;
    filtr_fd1(1, 2, 1, 1)=-1;

    //dmz/dz
    filtr_fd1(1, 1, 0, 2)= 1;
    filtr_fd1(1, 1, 2, 2)=-1;
    af::array first = convolve(state.m, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    //Make Divergence
    first=sum(first, 3);
    first=tile(first, 1, 1, 1, 3);
    first=n*first;

    //Dot product
    af::array second = sum(n*state.m, 3);
    //Expand for fd1 convolution
    second = tile(second, 1, 1, 1, 3);
    second = convolve(second, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    if(state.afsync){
        af::sync();
    }
    cpu_time += af::timer::stop(timer_dmi);
    return D_atom/(constants::mu0*state.Ms) * (first-second);
}

}// namespace magnumafcpp
