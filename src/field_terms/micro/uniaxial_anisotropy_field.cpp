// Ref Diss Abert p 15 sec 2.3 eq 2.3.28
// Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m)
// With higher order (not implemented): Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m) ( + 4
// Ku2 /(mu0 Ms) eu ( eu . m)^3
#include "micro/uniaxial_anisotropy_field.hpp"
#include "util/color_string.hpp" // color_string::warning()
#include "util/util.hpp"

namespace magnumafcpp {

UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, std::array<double, 3> Ku1_axis)
    : Ku1(Ku1), Ku1_axis(util::normalize_vector(Ku1_axis)) {}

UniaxialAnisotropyField::UniaxialAnisotropyField(af::array Ku1_field, std::array<double, 3> Ku1_axis)
    : Ku1_field(Ku1_field.dims(3) == 1 ? af::tile(Ku1_field, 1, 1, 1, 3) : std::move(Ku1_field)),
      Ku1_axis(util::normalize_vector(Ku1_axis)) {
    if (this->Ku1_field.dims(3) == 3) {
        printf("%s UniaxialAnisotropyField: You are using legacy dimension "
               "[nx, ny, nz, 3] for Ku1, please now use scalar field "
               "dimensions [nx, ny, nz, 1].\n",
               color_string::warning());
    }
}

// // Woulde be ambigous due to af::array non-explicit ctor
// UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, af::array Ku1_axis_field)
//     : Ku1(Ku1), Ku1_axis_field(util::normalize_handle_zero_vectors(Ku1_axis_field)) {}

UniaxialAnisotropyField::UniaxialAnisotropyField(af::array Ku1_field, const af::array& Ku1_axis_field)
    : Ku1_field(Ku1_field.dims(3) == 1 ? af::tile(Ku1_field, 1, 1, 1, 3) : std::move(Ku1_field)),
      Ku1_axis_field(util::normalize_handle_zero_vectors(Ku1_axis_field)) {
    if (this->Ku1_field.dims(3) == 3) {
        printf("%s UniaxialAnisotropyField: You are using legacy dimension "
               "[nx, ny, nz, 3] for Ku1, please now use scalar field "
               "dimensions [nx, ny, nz, 1].\n",
               color_string::warning());
    }
}

// // TODO Ku1 implicitly converted to af::array (as af::array ctor not explicit)!!!
// UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, long int Ku1_axis_field_ptr)
//     : UniaxialAnisotropyField(Ku1, util::pywrap::make_copy_form_py(Ku1_axis_field_ptr)) {
// }

// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, long int Ku1_axis_field_ptr)
    : Ku1(Ku1), Ku1_axis_field(util::normalize_handle_zero_vectors(util::pywrap::make_copy_form_py(Ku1_axis_field_ptr))) {}

// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField(long int Ku1_field_ptr, long int Ku1_axis_field_ptr)
    : UniaxialAnisotropyField(util::pywrap::make_copy_form_py(Ku1_field_ptr),
                              util::pywrap::make_copy_form_py(Ku1_axis_field_ptr)) {}

// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField(double Ku1, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2)
    : UniaxialAnisotropyField(Ku1, std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2}) {}

// For wrapping only
UniaxialAnisotropyField::UniaxialAnisotropyField(long int Ku1_field_ptr, double Ku1_axis_0, double Ku1_axis_1,
                                                 double Ku1_axis_2)
    : UniaxialAnisotropyField(util::pywrap::make_copy_form_py(Ku1_field_ptr),
                              std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2}) {}

af::array UniaxialAnisotropyField::impl_H_in_Apm(const State& state) const {
    // switch Ku1_axis and Ku1_axis_field
    af::array eu; // Array containing normal vectors
    if (Ku1_axis_field.isempty()) {
        eu = af::array(mesh::dims_v(state.mesh), f64);
        eu(af::span, af::span, af::span, 0) = Ku1_axis[0];
        eu(af::span, af::span, af::span, 1) = Ku1_axis[1];
        eu(af::span, af::span, af::span, 2) = Ku1_axis[2];
    } else {
        eu = Ku1_axis_field;
    }

    af::array anisotropy = eu * state.m;
    anisotropy = af::sum(anisotropy, 3);
    anisotropy = af::tile(anisotropy, 1, 1, 1, 3);

    if (state.Ms_field.isempty() && Ku1_field.isempty()) {
        return 2. * Ku1 / (constants::mu0 * state.Ms) * (eu * anisotropy);
    } else if (!state.Ms_field.isempty() && Ku1_field.isempty()) {
        af::array result = 2. * Ku1 / (constants::mu0 * state.Ms_field) * (eu * anisotropy);
        af::replace(result, !af::iszero(state.Ms_field),
                    0); // Replacing all resulting NaN with 0
        return result;
    } else if (state.Ms_field.isempty() && !Ku1_field.isempty()) {
        return 2. * Ku1_field / (constants::mu0 * state.Ms) * (eu * anisotropy);
    } else {
        af::array result = 2. * Ku1_field / (constants::mu0 * state.Ms_field) * (eu * anisotropy);
        af::replace(result, !af::iszero(state.Ms_field),
                    0); // Replacing all resulting NaN with 0
        return result;
    }
}

double UniaxialAnisotropyField::get_ku1_axis(int i) { return Ku1_axis[i]; }

long int UniaxialAnisotropyField::get_Ku1_field() const { return util::pywrap::send_copy_to_py(Ku1_field); }

} // namespace magnumafcpp
