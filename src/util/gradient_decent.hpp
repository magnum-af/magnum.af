#pragma once
#include "named_type.hpp"
#include<functional>
namespace magnumaf{

using X_start_val = NamedType<float, struct NamedTypeX_start_val>;
using Gamma = NamedType<float, struct NamedTypeGamma>;
using maxIters = NamedType<float, struct NamedTypemaxIters>;
using Epsilon = NamedType<float, struct NamedTypeEpsilon>;
using Precision = NamedType<float, struct NamedTypePrecision>;
using Verbose = NamedType<bool, struct NamedTypeVerbose>;
class GradientDecent{
    public:
        GradientDecent(std::function<float (float)> f, X_start_val x_start_val, Precision precision = Precision(1e-6), Gamma gamma = Gamma(1e-2), maxIters maxiters = maxIters(1000), Epsilon epsilon = Epsilon(1e-6), Verbose verbose = Verbose(true)):
            f(f), x_start_val(x_start_val), precision(precision), gamma(gamma), maxiters(maxiters), epsilon(epsilon), verbose(verbose)
    {
        //std::cout << precision.get() << std::endl;
        //std::cout << epsilon.get() << std::endl;
    }
        std::pair<float, float> minimize();
    private:
        std::function<float (float)> f;
        const X_start_val x_start_val;
        const Precision precision;
        const Gamma gamma;
        const maxIters maxiters;
        const Epsilon epsilon;
        float calc_df(float x_current);
        Verbose verbose;
};

}// namespace magnumaf
