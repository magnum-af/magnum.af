#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <utility>

namespace magnumafcpp {

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

} // namespace magnumafcpp
