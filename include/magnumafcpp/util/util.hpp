#pragma once
#include "arrayfire.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <utility>

namespace magnumafcpp {

namespace util::pywrap {

/// Copies an af::array from python.
/// Make an af::array from a pointer obtained from python
inline af::array make_copy_form_py(long int py_ptr_to_afarray) {
    return *(new af::array(*((void**)py_ptr_to_afarray)));
}

/// Passes a copy of the passed array to python.
/// Allocates a copy of passed af::array, returns the af_array handle as a long int
// alternative technical name: get_ptr_to_copy
inline long int send_copy_to_py(const af::array& afarray) { return (long int)(new af::array(afarray))->get(); }

} // namespace util::pywrap

/// template function calculating the mean and the standard deviation of a given
/// template container data standard deviation with -1
// from
// https://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
template <class T> std::pair<double, double> mean_stdev_w_minus(T vec) {
    double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);
    double m = sum / vec.size();
    double accum = 0.0;
    std::for_each(std::begin(vec), std::end(vec), [&](const double d) { accum += (d - m) * (d - m); });
    double stdev = std::sqrt(accum / (vec.size() - 1));
    return {m, stdev};
}

/// standard deviation without -1
template <class T> std::pair<double, double> mean_stdev_no_minus(T vec) {
    double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);
    double m = sum / vec.size();
    double accum = 0.0;
    std::for_each(std::begin(vec), std::end(vec), [&](const double d) { accum += (d - m) * (d - m); });
    double stdev = std::sqrt(accum / (vec.size()));
    return {m, stdev};
}

inline std::array<double, 3> cross_product(std::array<double, 3> a, std::array<double, 3> b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

inline double dot_product(std::array<double, 3> a, std::array<double, 3> b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline double vector_norm(std::array<double, 3> vector) {
    return std::sqrt(std::pow(vector[0], 2) + std::pow(vector[1], 2) + std::pow(vector[2], 2));
}

inline std::array<double, 3> normalize_vector(std::array<double, 3> vector) {
    double norm = vector_norm(vector);
    return {vector[0] / norm, vector[1] / norm, vector[2] / norm};
}

} // namespace magnumafcpp
