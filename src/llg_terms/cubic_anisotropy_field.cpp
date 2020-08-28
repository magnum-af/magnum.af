#include "cubic_anisotropy_field.hpp"
#include "util.hpp"
namespace magnumafcpp {

CubicAnisotropyField::CubicAnisotropyField(double Kc1, double Kc2, double Kc3, std::array<double, 3> c1,
                                           std::array<double, 3> c2)
    : Kc1(Kc1), Kc2(Kc2), Kc3(Kc3), c1(normalize_vector(c1)), c2(normalize_vector(c2)),
      c3(normalize_vector(cross_product(c1, c2))) {
    if (dot_product(c1, c2) > 0) {
        std::cout << "Warning in CubicAnisotropyField: provided c1 and c2 are not perpendicular, (c1 . c2) > 0"
                  << std::endl;
        std::cout << "Please choose perpendicular input vectors." << std::endl;
        exit(1);
    }
}

CubicAnisotropyField::CubicAnisotropyField(double Kc1, double Kc2, double Kc3, double c1x, double c1y, double c1z,
                                           double c2x, double c2y, double c2z)
    : CubicAnisotropyField(Kc1, Kc2, Kc3, {c1x, c1y, c1z}, {c2x, c2y, c2z}) {}

af::array dot_4d(const af::array& a, const af::array& b) { return af::tile(af::sum(a * b, 3), 1, 1, 1, 3); }

std::array<af::array, 3> CubicAnisotropyField::h_1to3(const State& state) {
    // c1,c2,c3 are double so c1_, c2_, c3_ initially are f64 before .as()
    af::array c1_ = af::tile(af::array(1, 1, 1, 3, c1.data()).as(state.m.type()), state.m.dims(0), state.m.dims(1),
                             state.m.dims(2), 1);
    af::array c2_ = af::tile(af::array(1, 1, 1, 3, c2.data()).as(state.m.type()), state.m.dims(0), state.m.dims(1),
                             state.m.dims(2), 1);
    af::array c3_ = af::tile(af::array(1, 1, 1, 3, c3.data()).as(state.m.type()), state.m.dims(0), state.m.dims(1),
                             state.m.dims(2), 1);

    af::array c1m = dot_4d(c1_, state.m);
    af::array c2m = dot_4d(c2_, state.m);
    af::array c3m = dot_4d(c3_, state.m);
    af::array c1m2 = af::pow(c1m, 2);
    af::array c2m2 = af::pow(c2m, 2);
    af::array c3m2 = af::pow(c3m, 2);

    af::array h1 = -2 * Kc1 / (constants::mu0 * state.Ms) *
                   ((c2m2 + c3m2) * c1m * c1_ + (c1m2 + c3m2) * c2m * c2_ + (c1m2 + c2m2) * c3m * c3_);

    af::array h2 = -2 * Kc2 / (constants::mu0 * state.Ms) *
                   (c2m2 * c3m2 * c1m * c1_ + c1m2 * c3m2 * c2m * c2_ + c1m2 * c2m2 * c3m * c3_);

    af::array c1m4 = af::pow(c1m, 4);
    af::array c2m4 = af::pow(c2m, 4);
    af::array c3m4 = af::pow(c3m, 4);
    af::array c1m3 = af::pow(c1m, 3);
    af::array c2m3 = af::pow(c2m, 3);
    af::array c3m3 = af::pow(c3m, 3);
    af::array h3 = -4 * Kc3 / (constants::mu0 * state.Ms) *
                   ((c2m4 + c3m4) * c1m3 * c1_ + (c1m4 + c3m4) * c2m3 * c2_ + (c1m4 + c2m4) * c3m3 * c3_);
    return {h1, h2, h3};
}

af::array CubicAnisotropyField::h(const State& state) {
    auto h = h_1to3(state);
    return h[0] + h[1] + h[2];
}

double CubicAnisotropyField::E(const State& state) {
    auto h = h_1to3(state);
    return constants::mu0 * state.Ms *
           af::sum(af::sum(af::sum(af::sum((-1 / 4. * h[0] - 1 / 6. * h[1] - 1 / 8. * h[2]) * state.m, 0), 1), 2), 3)
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
