#include "micro/spin_transfer_torque_field.hpp"
#include "math.hpp"

namespace magnumafcpp {

SpinTransferTorqueField::SpinTransferTorqueField(af::array polarization_field, double nu_dampinglike,
                                                 double nu_fieldlike, double j_e)
    : polarization_field(polarization_field), nu_dampinglike(nu_dampinglike), nu_fieldlike(nu_fieldlike), j_e(j_e) {}

SpinTransferTorqueField::SpinTransferTorqueField(long int polarization_field_ptr, double nu_dampinglike,
                                                 double nu_fieldlike, double j_e)
    : polarization_field(polarization_field_ptr), nu_dampinglike(nu_dampinglike), nu_fieldlike(nu_fieldlike), j_e(j_e) {
}

af::array SpinTransferTorqueField::impl_H_in_Apm(const State& state) const {
    if (state.Ms_field.isempty()) {
        return -j_e * constants::hbar / (2. * constants::e_neg * constants::mu0 * state.Ms) *
               (nu_dampinglike * math::cross4(state.m, polarization_field.array) +
                nu_fieldlike * polarization_field.array);
    } else {
        af::array h_sot = -j_e * constants::hbar / (2. * constants::e_neg * constants::mu0 * state.Ms_field) *
                          (nu_dampinglike * math::cross4(state.m, polarization_field.array) +
                           nu_fieldlike * polarization_field.array);
        af::replace(h_sot, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return h_sot;
    }
}
} // namespace magnumafcpp
