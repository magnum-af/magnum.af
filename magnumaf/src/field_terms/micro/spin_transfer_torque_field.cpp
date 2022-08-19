#include "micro/spin_transfer_torque_field.hpp"
#include "math.hpp"
#include "util/util.hpp"

namespace magnumaf {

SpinTransferTorqueField::SpinTransferTorqueField(const af::array& polarization, double eta_damping, double eta_field,
                                                 double j_e, double fl_thickness)
    : polarization(util::normalize(polarization)), eta_damping(eta_damping), eta_field(eta_field), j_e(j_e),
      fl_thickness(fl_thickness) {}

SpinTransferTorqueField::SpinTransferTorqueField(long int polarization_ptr, double eta_damping, double eta_field,
                                                 double j_e, double fl_thickness)
    : SpinTransferTorqueField(util::pywrap::make_copy_form_py(polarization_ptr), eta_damping, eta_field, j_e,
                              fl_thickness) {}

af::array SpinTransferTorqueField::impl_H_in_Apm(const State& state) const {
    if (state.Ms_field.isempty()) {
        return -j_e * constants::hbar / (2. * constants::e_neg * constants::mu0 * state.Ms * fl_thickness) *
               (eta_damping * math::cross4(state.m, polarization) + eta_field * polarization);
    } else {
        af::array h_sot = -j_e * constants::hbar /
                          (2. * constants::e_neg * constants::mu0 * state.get_Ms_field_in_vec_dims() * fl_thickness) *
                          (eta_damping * math::cross4(state.m, polarization) + eta_field * polarization);
        af::replace(h_sot, state.get_Ms_field_in_vec_dims() != 0, 0); // set all cells where Ms==0 to 0
        return h_sot;
    }
}
} // namespace magnumaf
