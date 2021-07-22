#include "nonequi/nonequi_uniaxial_anisotropy_field.hpp"


#include <utility>

#include "util/util.hpp"

namespace magnumaf {

NonequiUniaxialAnisotropyField::NonequiUniaxialAnisotropyField(NonequiMesh nemesh, double Ku1,
                                                               std::array<double, 3> Ku1_axis)
    : NonequiTerm(std::move(nemesh)), Ku1(Ku1), Ku1_axis(get_normalized_vector(Ku1_axis)) {}

NonequiUniaxialAnisotropyField::NonequiUniaxialAnisotropyField(NonequiMesh nemesh, af::array Ku1_field,
                                                               std::array<double, 3> Ku1_axis)
    : NonequiTerm(std::move(nemesh)), Ku1_field(std::move(Ku1_field)), Ku1_axis(get_normalized_vector(Ku1_axis)) {}

NonequiUniaxialAnisotropyField::NonequiUniaxialAnisotropyField(NonequiMesh nemesh, af::array Ku1_field,
                                                               af::array Ku1_axis_field)
    : NonequiTerm(std::move(nemesh)), Ku1_field(std::move(Ku1_field)), Ku1_axis_field(std::move(Ku1_axis_field)) {}

NonequiUniaxialAnisotropyField::NonequiUniaxialAnisotropyField(NonequiMesh nemesh, double Ku1,
                                                               long int Ku1_axis_field_ptr)
    : NonequiTerm(std::move(nemesh)), Ku1(Ku1), Ku1_axis_field(util::pywrap::make_copy_form_py(Ku1_axis_field_ptr)) {}

NonequiUniaxialAnisotropyField::NonequiUniaxialAnisotropyField(NonequiMesh nemesh, long int Ku1_field_ptr,
                                                               long int Ku1_axis_field_ptr)
    : NonequiTerm(std::move(nemesh)), Ku1_field(util::pywrap::make_copy_form_py(Ku1_field_ptr)),
      Ku1_axis_field(util::pywrap::make_copy_form_py(Ku1_axis_field_ptr)) {}

// For wrapping only
NonequiUniaxialAnisotropyField::NonequiUniaxialAnisotropyField(NonequiMesh nemesh, double Ku1, double Ku1_axis_0,
                                                               double Ku1_axis_1, double Ku1_axis_2)
    : NonequiTerm(std::move(nemesh)), Ku1(Ku1),
      Ku1_axis(get_normalized_vector(std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2})) {}

// For wrapping only
NonequiUniaxialAnisotropyField::NonequiUniaxialAnisotropyField(NonequiMesh nemesh, long int Ku1_field_ptr,
                                                               double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2)
    : NonequiTerm(std::move(nemesh)), Ku1_field(util::pywrap::make_copy_form_py(Ku1_field_ptr)),
      Ku1_axis(get_normalized_vector(std::array<double, 3>{Ku1_axis_0, Ku1_axis_1, Ku1_axis_2})) {}

af::array NonequiUniaxialAnisotropyField::impl_H_in_Apm(const State& state) const { return calc_heff(state); }

af::array NonequiUniaxialAnisotropyField::calc_heff(const State& state) const {
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

std::array<double, 3> NonequiUniaxialAnisotropyField::get_normalized_vector(std::array<double, 3> vector) {
    double norm = sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
    return std::array<double, 3>{vector[0] / norm, vector[1] / norm, vector[2] / norm};
}

double NonequiUniaxialAnisotropyField::get_ku1_axis(int i) const { return Ku1_axis[i]; }

long int NonequiUniaxialAnisotropyField::get_Ku1_field() const { return util::pywrap::send_copy_to_py(Ku1_field); }
} // namespace magnumaf
