#include "math.hpp"
namespace magnumafcpp::math{

double max_4d_abs(const af::array& a) {
    return af::max(af::max(af::max(af::max(af::abs(a), 0), 1), 2), 3).as(f64).scalar<double>();
}

af::array cross4(const af::array& a, const af::array& b) {
    af::array c = af::array(a.dims(0), a.dims(1), a.dims(2), 3, a.type());
    c(af::span, af::span, af::span, 0) = a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 2) -
                                         a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 1);
    c(af::span, af::span, af::span, 1) = a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 0) -
                                         a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 2);
    c(af::span, af::span, af::span, 2) = a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 1) -
                                         a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 0);
    return c;
}

af::array cross4shift(const af::array& a, const af::array& b) {
    af::array ashift = af::shift(a, 0, 0, 0, -1);
    af::array ashift2 = af::shift(a, 0, 0, 0, -2);
    af::array bshift = af::shift(b, 0, 0, 0, -2);
    af::array bshift2 = af::shift(b, 0, 0, 0, -1);
    return ashift * bshift - ashift2 * bshift2;
}

} // namespace magnumafcpp::math
