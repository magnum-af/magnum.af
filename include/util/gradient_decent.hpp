#pragma once
#include "named_type.hpp"
#include <functional>
namespace magnumafcpp
{

using X_start_val = NamedType<double, struct NamedTypeX_start_val>;
using Gamma = NamedType<double, struct NamedTypeGamma>;
using maxIters = NamedType<double, struct NamedTypemaxIters>;
using Epsilon = NamedType<double, struct NamedTypeEpsilon>;
using Precision = NamedType<double, struct NamedTypePrecision>;
using Verbose = NamedType<bool, struct NamedTypeVerbose>;
class GradientDecent
{
public:
    GradientDecent(std::function<double(double)> f, X_start_val x_start_val, Precision precision = Precision(1e-6), Gamma gamma = Gamma(1e-2), maxIters maxiters = maxIters(1000), Epsilon epsilon = Epsilon(1e-6), Verbose verbose = Verbose(true)) : f(f), x_start_val(x_start_val), precision(precision), gamma(gamma), maxiters(maxiters), epsilon(epsilon), verbose(verbose)
    {
        //std::cout << precision.get() << std::endl;
        //std::cout << epsilon.get() << std::endl;
    }
    std::pair<double, double> minimize();

private:
    std::function<double(double)> f;
    const X_start_val x_start_val;
    const Precision precision;
    const Gamma gamma;
    const maxIters maxiters;
    const Epsilon epsilon;
    double calc_df(double x_current);
    Verbose verbose;
};

} // namespace magnumafcpp
