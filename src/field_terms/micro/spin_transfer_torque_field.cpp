#include "micro/spin_transfer_torque_field.hpp"
#include "math.hpp"
#include "util/util.hpp"

namespace magnumafcpp {

SpinTransferTorqueField::SpinTransferTorqueField(af::array polarization, double nu_damping, double nu_field, double j_e,
                                                 double fl_thickness)
    : polarization(std::move(polarization)), nu_damping(nu_damping), nu_field(nu_field), j_e(j_e),
      fl_thickness(fl_thickness) {}

SpinTransferTorqueField::SpinTransferTorqueField(long int polarization_ptr, double nu_damping, double nu_field,
                                                 double j_e, double fl_thickness)
    : SpinTransferTorqueField(util::pywrap::make_copy_form_py(polarization_ptr), nu_damping, nu_field, j_e,
                              fl_thickness) {}

af::array SpinTransferTorqueField::impl_H_in_Apm(const State& state) const {
    if (state.Ms_field.isempty()) {
        return -j_e * constants::hbar / (2. * constants::e_neg * constants::mu0 * state.Ms * fl_thickness) *
               (nu_damping * math::cross4(state.m, polarization) + nu_field * polarization);
    } else {
        af::array h_sot = -j_e * constants::hbar /
                          (2. * constants::e_neg * constants::mu0 * state.Ms_field * fl_thickness) *
                          (nu_damping * math::cross4(state.m, polarization) + nu_field * polarization);
        af::replace(h_sot, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return h_sot;
    }
}
} // namespace magnumafcpp
