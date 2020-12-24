#pragma once
#include<array>
#include"arrayfire.h"
namespace magnumafcpp::math{

/// Calculate mean along first three dimensions of af::array.
template <typename T> std::array<T, 3> mean_3d(const af::array& a) {
    const auto mean = af::mean(af::mean(af::mean(a, 0), 1), 2);
    const double mx = mean(0, 0, 0, 0).scalar<T>();
    const double my = mean(0, 0, 0, 1).scalar<T>();
    const double mz = mean(0, 0, 0, 2).scalar<T>();
    return {mx, my, mz};
}

/// Absolute value of maximum of all values in array
inline double max_4d_abs(const af::array& a) {
    return af::max(af::max(af::max(af::max(af::abs(a), 0), 1), 2), 3).as(f64).scalar<double>();
}

/// Cross product for vector fields of format [nx, ny, nz, 3]
inline af::array cross4(const af::array& a, const af::array& b) {
    af::array c = af::array(a.dims(0), a.dims(1), a.dims(2), 3, a.type());
    c(af::span, af::span, af::span, 0) = a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 2) -
                                         a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 1);
    c(af::span, af::span, af::span, 1) = a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 0) -
                                         a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 2);
    c(af::span, af::span, af::span, 2) = a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 1) -
                                         a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 0);
    return c;
}

// Slightly slower than cross4
/// Cross product for vector fields of format [nx, ny, nz, 3] using af::shift
inline af::array cross4shift(const af::array& a, const af::array& b) {
    af::array ashift = af::shift(a, 0, 0, 0, -1);
    af::array ashift2 = af::shift(a, 0, 0, 0, -2);
    af::array bshift = af::shift(b, 0, 0, 0, -2);
    af::array bshift2 = af::shift(b, 0, 0, 0, -1);
    return ashift * bshift - ashift2 * bshift2;
}

} // namespace magnumafcpp::math
