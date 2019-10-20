#pragma once
#include "named_type.hpp"
#include<functional>
#include<math.h>
#include<iostream>
namespace magnumaf{

using X0 = NamedType<float, struct NamedTypeX0>;
using Imax = NamedType<unsigned, struct NamedTypeImax>;
using Precision = NamedType<float, struct NamedTypePrecision>;
using EpsilonFactor = NamedType<float, struct NamedTypeEpsilonFactor>;
using Verbose = NamedType<bool, struct NamedTypeVerbose>;

/// Netwon Iteration Object

/// Tries to find x s.t. f(x) = 0.
/// Returns x and f(x) as std::pair <float, float> {x, f(x)}.
/// \param EpsilonFactor used for adaptively calculating epsilon \f$ eps = EpsilonFactor * |x_n - x_{n-1}|\f$ for the approximation of the derivative: \f$ f'(x) = \frac{f(x+eps) - f(x)}{eps} \f$

class NewtonIteration{
    public:
        NewtonIteration(std::function<float (float)> f, Verbose verbose = Verbose(true)): f(f), verbose(verbose){}

        std::pair<float, float> run (X0 x, Precision precision = Precision(1e-8), EpsilonFactor epsfac = EpsilonFactor(1e-4), Imax imax = Imax(100)){
            float eps_factor_wrt_xrange = 1e-3;
            float eps = epsfac.get() * x.get();
            for (unsigned i = 0; i < imax.get(); i ++){
                auto f_and_df = df(x.get(), eps);
                float x_current = x.get();
                x.get() -= f_and_df.first/f_and_df.second;
                eps =  eps_factor_wrt_xrange * fabs(x.get() - x_current);
                if (verbose.get()) std::cout << "i=" << i << ", x=" << x.get() << ", xold=" << x_current <<", f(x)=" << f_and_df.first << ", df(x)=" << f_and_df.second << ", eps=" << eps << std::endl;
                // return if f(x) < precision
                if (std::fabs(f_and_df.first) < precision.get()) {
                    if (verbose.get()) printf("NewtonIteration.run: precision reached, finished.\n");
                    return std::pair<float, float> (x_current, f_and_df.first);
                }
                // abort if eps approaches 0 or inf
                if (eps < 1e-100 or eps > 1e100) {
                    printf("NewtonIteration: eps=%f is < 1e-100 or > 1e100, abort.\n", eps);
                    return std::pair<float, float> (x_current, f_and_df.first);
                }
            }
            // return if imax is reached
            return std::pair<float, float> (x.get(), f(x.get()));
        }
    private:
        // calculate f(x) and df(x) and return as std::pair <float, float> {f(x), df(x)}
        std::pair<float, float> df(float x, float eps){
                float f_x_eps = f(x + eps);
                float f_x = f(x);
                float df_x = ( f_x_eps - f_x )/eps;
                return {f_x, df_x};
        }

        std::function<float (float)> f;
        Verbose verbose;
};

}// namespace magnumaf
