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

} // namespace magnumafcpp::math
