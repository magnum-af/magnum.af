#include "cubic_anisotropy_field.hpp"
#include "func.hpp"
#include "util.hpp"
namespace magnumafcpp {

af::array dot_4d(const af::array& a, const af::array& b) { return af::tile(af::sum(a * b, 3), 1, 1, 1, 3); }

CubicAnisotropyField::CubicAnisotropyField(double Kc1, double Kc2, double Kc3, std::array<double, 3> c1,
                                           std::array<double, 3> c2)
    : Kc1(Kc1), Kc2(Kc2), Kc3(Kc3), c1(normalize_vector(c1)), c2(normalize_vector(c2)),
      c3(normalize_vector(cross_product(c1, c2))) {
    // check input vectors c1, c2
    const double precision = 1e-12;
    const double abs_dot_c1c2 = std::fabs(dot_product(this->c1, this->c2));
    if (abs_dot_c1c2 > precision) {
        std::cout << "Warning in CubicAnisotropyField: provided c1 and c2 are not perpendicular, i.e. (c1 . c2) = "
                  << abs_dot_c1c2 << " > 0" << std::endl;
        std::cout << "Please choose perpendicular input vectors." << std::endl;
        exit(1);
    }
}

CubicAnisotropyField::CubicAnisotropyField(af::array Kc1_array, af::array Kc2_array, af::array Kc3_array,
                                           af::array c1_array, af::array c2_array)
    : Kc1_array(Kc1_array), Kc2_array(Kc2_array), Kc3_array(Kc3_array),
      c1_array(normalize_handle_zero_vectors(c1_array)), c2_array(normalize_handle_zero_vectors(c2_array)),
      c3_array(cross4(this->c1_array, this->c2_array)) {
    // check input vectors c1, c2
    const double precision = 1e-12;
    const double max_abs_c1_c2_dot =
        af::max(af::max(af::max(af::max(af::abs(dot_4d(this->c1_array, this->c2_array)), 0), 1), 2), 3)
            .scalar<double>();
    if (max_abs_c1_c2_dot > precision) {
        std::cout << "Warning in CubicAnisotropyField: provided c1 and c2 are not perpendicular, i.e. "
                     "max(abs((c1 . c2)) = "
                  << max_abs_c1_c2_dot << " > " << precision << std::endl;
        std::cout << "Please choose perpendicular input vectors." << std::endl;
        exit(1);
    }
}

af::array ptr_to_array(long int array_ptr) { return *(new af::array(*((void**)array_ptr))); }

// Wrapping
CubicAnisotropyField::CubicAnisotropyField(long int Kc1_array_ptr, long int Kc2_array_ptr, long int Kc3_array_ptr,
                                           long int c1_array_ptr, long int c2_array_ptr)
    : CubicAnisotropyField(ptr_to_array(Kc1_array_ptr), ptr_to_array(Kc2_array_ptr), ptr_to_array(Kc3_array_ptr),
                           ptr_to_array(c1_array_ptr), ptr_to_array(c2_array_ptr)) {}

// Wrapping
CubicAnisotropyField::CubicAnisotropyField(double Kc1, double Kc2, double Kc3, double c1x, double c1y, double c1z,
                                           double c2x, double c2y, double c2z)
    : CubicAnisotropyField(Kc1, Kc2, Kc3, {c1x, c1y, c1z}, {c2x, c2y, c2z}) {}


std::array<af::array, 3> CubicAnisotropyField::h_1to3(const State& state) {
    // c1,c2,c3 are double so c1_, c2_, c3_ initially are f64 before .as()
    af::array Kc1_, Kc2_, Kc3_;

    if (Kc1_array.isempty()) { // Assuming Kc*_array are all empty

        Kc1_ = af::constant(Kc1, state.m.dims(), state.m.type());
        Kc2_ = af::constant(Kc2, state.m.dims(), state.m.type());
        Kc3_ = af::constant(Kc3, state.m.dims(), state.m.type());
    } else {
        Kc1_ = af::tile(Kc1_array, 1, 1, 1, 3);
        Kc2_ = af::tile(Kc2_array, 1, 1, 1, 3);
        Kc3_ = af::tile(Kc3_array, 1, 1, 1, 3);
    }

    af::array c1_, c2_, c3_;
    if (c1_array.isempty()) {
        c1_ = af::tile(af::array(1, 1, 1, 3, c1.data()).as(state.m.type()), state.m.dims(0), state.m.dims(1),
                       state.m.dims(2), 1);
        c2_ = af::tile(af::array(1, 1, 1, 3, c2.data()).as(state.m.type()), state.m.dims(0), state.m.dims(1),
                       state.m.dims(2), 1);
        c3_ = af::tile(af::array(1, 1, 1, 3, c3.data()).as(state.m.type()), state.m.dims(0), state.m.dims(1),
                       state.m.dims(2), 1);
    } else {
        c1_ = c1_array;
        c2_ = c2_array;
        c3_ = c3_array;
    }

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

    const af::array Ms_ = state.get_Ms_field_in_vector_dims();

    af::array h1 = -2 * Kc1_ / (constants::mu0 * Ms_) *
                   ((c2m2 + c3m2) * c1m * c1_ + (c1m2 + c3m2) * c2m * c2_ + (c1m2 + c2m2) * c3m * c3_);

    af::array h2 = -2 * Kc2_ / (constants::mu0 * Ms_) *
                   (c2m2 * c3m2 * c1m * c1_ + c1m2 * c3m2 * c2m * c2_ + c1m2 * c2m2 * c3m * c3_);

    af::array h3 = -4 * Kc3_ / (constants::mu0 * Ms_) *
                   ((c2m4 + c3m4) * c1m3 * c1_ + (c1m4 + c3m4) * c2m3 * c2_ + (c1m4 + c2m4) * c3m3 * c3_);
    return {h1, h2, h3};
}

af::array CubicAnisotropyField::h(const State& state) {
    auto h = h_1to3(state);
    return h[0] + h[1] + h[2];
}

double CubicAnisotropyField::E(const State& state) {
    const af::array Ms = state.get_Ms_field_in_vector_dims();
    auto h = h_1to3(state);
    return constants::mu0 *
           af::sum(af::sum(af::sum(af::sum(Ms * (-1 / 4. * h[0] - 1 / 6. * h[1] - 1 / 8. * h[2]) * state.m, 0), 1), 2),
                   3)
               .as(f64)
               .scalar<double>() *
           state.mesh.dx * state.mesh.dy * state.mesh.dz;
}

// TODO reglects h caching:
double CubicAnisotropyField::E(const State& state, const af::array& h) {
    auto h_avoid_warning = h;
    return E(state);
}

} // namespace magnumafcpp
