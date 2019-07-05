#pragma once
#include <numeric>
#include <utility>
#include <algorithm>

namespace magnumaf{

/// template function calculating the mean and standard deviation of a given container data
// from https://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
template<class T>
std::pair<double, double> mean_stdev(T vec){
    double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);
    double m =  sum / vec.size();
    double accum = 0.0;
    std::for_each (std::begin(vec), std::end(vec), [&](const double d) {
        accum += (d - m) * (d - m);
    });
    double stdev = sqrt(accum / (vec.size()-1));
    return {m, stdev};
}

}// namespace magnumaf
