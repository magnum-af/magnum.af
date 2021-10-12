#include "field_terms/micro/dmi_d2d_field.hpp"
namespace magnumaf {

af::array DMI_D2d_Field::impl_H_in_Apm(const State& state) const {
    // neighbour in x+/x-/y+/y-:
    af::array x_plus = af::shift(state.m, 1) / state.mesh.dx;
    af::array x_minus = af::shift(state.m, -1) / state.mesh.dx;
    af::array y_plus = af::shift(state.m, 0, 1) / state.mesh.dy;
    af::array y_minus = af::shift(state.m, 0, -1) / state.mesh.dy;

    // setting non-existant neighbours to zero if not in PBC mode:
    if (not PBC_) {
        x_plus(0, af::span, af::span, af::span) = 0.0;
        x_minus(-1, af::span, af::span, af::span) = 0.0;
        y_plus(af::span, 0, af::span, af::span) = 0.0;
        y_minus(af::span, -1, af::span, af::span) = 0.0;
    }

    // selecting mx/my/mz components:
    const af::array y_plus_mx = y_plus(af::span, af::span, af::span, 0);
    const af::array y_plus_mz = y_plus(af::span, af::span, af::span, 2);

    const af::array y_minus_mx = y_minus(af::span, af::span, af::span, 0);
    const af::array y_minus_mz = y_minus(af::span, af::span, af::span, 2);

    const af::array x_minus_mz = x_minus(af::span, af::span, af::span, 2);
    const af::array x_minus_my = x_minus(af::span, af::span, af::span, 1);

    const af::array x_plus_mz = x_plus(af::span, af::span, af::span, 2);
    const af::array x_plus_my = x_plus(af::span, af::span, af::span, 1);

    // composing result of cross-product with unit vectors:
    const af::array Heff_x = y_plus_mz - y_minus_mz;
    const af::array Heff_y = x_plus_mz - x_minus_mz;
    const af::array Heff_z = y_minus_mx - y_plus_mx + x_minus_my - x_plus_my;
    af::array Heff = af::join(3, Heff_x, Heff_y, Heff_z);

    // Adding global coefficients:
    Heff *= -D_in_J_per_m2_ / (constants::mu0);
    if (state.Ms_field.isempty()) {
        return Heff / state.Ms;
    } else {
        af::array result = Heff / state.Ms_field;
        af::replace(result, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return result;
    }
}
} // namespace magnumaf
