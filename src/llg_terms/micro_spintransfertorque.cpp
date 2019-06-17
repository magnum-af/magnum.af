#include "micro_spintransfertorque.hpp"

double SpinTransferTorqueField::E(const State& state){
    return - constants::mu0/2. * state.Ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}

double SpinTransferTorqueField::E(const State& state, const af::array& h){
    return -constants::mu0/2. * state.Ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}


SpinTransferTorqueField::SpinTransferTorqueField (af::array polarization_field, double nu_dampinglike, double nu_fieldlike, double j_e) : polarization_field(polarization_field), nu_dampinglike(nu_dampinglike), nu_fieldlike(nu_fieldlike), j_e(j_e) {
}

SpinTransferTorqueField::SpinTransferTorqueField (long int polarization_field_ptr, double nu_dampinglike, double nu_fieldlike, double j_e) : polarization_field(polarization_field_ptr), nu_dampinglike(nu_dampinglike), nu_fieldlike(nu_fieldlike), j_e(j_e) {
}

af::array SpinTransferTorqueField::h(const State& state){
    //af::timer timer_fieldlike = af::timer::start();
    //evaluation_timing += af::timer::stop(timer_fieldlike);
    return - j_e * constants::hbar / (2. * constants::e * constants::mu0 * state.Ms) * (nu_dampinglike * cross4(state.m, polarization_field.array) + nu_fieldlike * polarization_field.array);
}
