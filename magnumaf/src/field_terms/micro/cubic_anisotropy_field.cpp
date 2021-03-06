#include "field_terms/micro/cubic_anisotropy_field.hpp"
#include "math.hpp"
#include "util/util.hpp"
#include <variant>
namespace magnumaf {

af::array dot_4d(const af::array& a, const af::array& b) { return af::tile(af::sum(a * b, 3), 1, 1, 1, 3); }

CubicAnisotropyField::CubicAnisotropyField(double Kc1, double Kc2, double Kc3, std::array<double, 3> c1,
                                           std::array<double, 3> c2)
    : Kc1(Kc1), Kc2(Kc2), Kc3(Kc3), c1(c1), c2(c2), c3(util::normalize_vector(util::cross_product(c1, c2))) {
    // check input vectors c1, c2
    const double precision = 1e-12;
    const double abs_dot_c1c2 = std::fabs(util::dot_product(std::get<std::array<double, 3>>(this->c1.variant),
                                                            std::get<std::array<double, 3>>(this->c2.variant)));
    if (abs_dot_c1c2 > precision) {
        std::cout << "Warning in CubicAnisotropyField: provided c1 and c2 are not perpendicular, i.e. (c1 . c2) = "
                  << abs_dot_c1c2 << " > 0" << std::endl;
        std::cout << "Please choose perpendicular input vectors." << std::endl;
        exit(1);
    }
}

CubicAnisotropyField::CubicAnisotropyField(af::array Kc1_array, af::array Kc2_array, af::array Kc3_array, af::array c1,
                                           af::array c2)
    : Kc1(std::move(Kc1_array)), Kc2(std::move(Kc2_array)), Kc3(std::move(Kc3_array)), c1(std::move(c1)),
      c2(std::move(c2)),
      c3(math::cross4(std::get<af::array>(this->c1.variant), std::get<af::array>(this->c2.variant))) {
    // check input vectors c1, c2
    const double precision = 1e-12;
    const auto max_abs_c1_c2_dot =
        af::max(af::max(af::max(af::max(af::abs(dot_4d(std::get<af::array>(this->c1.variant),
                                                       std::get<af::array>(this->c2.variant))),
                                        0),
                                1),
                        2),
                3)
            .scalar<double>();
    if (max_abs_c1_c2_dot > precision) {
        std::cout << "Warning in CubicAnisotropyField: provided c1 and c2 are not perpendicular, i.e. "
                     "max(abs((c1 . c2)) = "
                  << max_abs_c1_c2_dot << " > " << precision << std::endl;
        std::cout << "Please choose perpendicular input vectors." << std::endl;
        exit(1);
    }
}

af::array ptr_to_array(long int array_ptr) { return util::pywrap::make_copy_form_py(array_ptr); }

// Wrapping
CubicAnisotropyField::CubicAnisotropyField(long int Kc1_array_ptr, long int Kc2_array_ptr, long int Kc3_array_ptr,
                                           long int c1_array_ptr, long int c2_array_ptr)
    : CubicAnisotropyField(ptr_to_array(Kc1_array_ptr), ptr_to_array(Kc2_array_ptr), ptr_to_array(Kc3_array_ptr),
                           ptr_to_array(c1_array_ptr), ptr_to_array(c2_array_ptr)) {}

// Wrapping
CubicAnisotropyField::CubicAnisotropyField(double Kc1, double Kc2, double Kc3, double c1x, double c1y, double c1z,
                                           double c2x, double c2y, double c2z)
    : CubicAnisotropyField(Kc1, Kc2, Kc3, {c1x, c1y, c1z}, {c2x, c2y, c2z}) {}

std::array<af::array, 3> CubicAnisotropyField::h_1to3(const State& state) const {
    af::array c1_ = c1.get_as_array(state.m.dims(), state.m.type());
    af::array c2_ = c2.get_as_array(state.m.dims(), state.m.type());
    af::array c3_ = c3.get_as_array(state.m.dims(), state.m.type());

    af::array c1m = dot_4d(c1_, state.m);
    af::array c2m = dot_4d(c2_, state.m);
    af::array c3m = dot_4d(c3_, state.m);
    af::array c1m2 = af::pow(c1m, 2);
    af::array c2m2 = af::pow(c2m, 2);
    af::array c3m2 = af::pow(c3m, 2);

    af::array c1m4 = af::pow(c1m, 4);
    af::array c2m4 = af::pow(c2m, 4);
    af::array c3m4 = af::pow(c3m, 4);
    af::array c1m3 = af::pow(c1m, 3);
    af::array c2m3 = af::pow(c2m, 3);
    af::array c3m3 = af::pow(c3m, 3);

    const af::array Ms_ = state.get_Ms_as_field_in_vector_dims();

    af::array h1 = -2 / (constants::mu0 * Ms_) * Kc1 *
                   ((c2m2 + c3m2) * c1m * c1_ + (c1m2 + c3m2) * c2m * c2_ + (c1m2 + c2m2) * c3m * c3_);

    af::array h2 = -2 / (constants::mu0 * Ms_) * Kc2 *
                   (c2m2 * c3m2 * c1m * c1_ + c1m2 * c3m2 * c2m * c2_ + c1m2 * c2m2 * c3m * c3_);

    af::array h3 = -4 / (constants::mu0 * Ms_) * Kc3 *
                   ((c2m4 + c3m4) * c1m3 * c1_ + (c1m4 + c3m4) * c2m3 * c2_ + (c1m4 + c2m4) * c3m3 * c3_);
    return {h1, h2, h3};
}

af::array CubicAnisotropyField::impl_H_in_Apm(const State& state) const {
    auto h = h_1to3(state);
    return h[0] + h[1] + h[2];
}

double CubicAnisotropyField::impl_E_in_J(const State& state, const af::array& h) const {
    // Note, h is ignored here, we need h_1to3
    // would require interface exception
    h.isempty(); // avoiding unused warning
    const af::array Ms = state.get_Ms_as_field_in_vector_dims();
    auto htemp = h_1to3(state);
    return constants::mu0 *
           af::sum(
               af::sum(
                   af::sum(af::sum(Ms * (-1 / 4. * htemp[0] - 1 / 6. * htemp[1] - 1 / 8. * htemp[2]) * state.m, 0), 1),
                   2),
               3)
               .as(f64)
               .scalar<double>() *
           state.mesh.dx * state.mesh.dy * state.mesh.dz;
}

} // namespace magnumaf
