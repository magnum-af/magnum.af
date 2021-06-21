#include "bs23.hpp"
#include "arrayfire.h"
#include "math.h"

#include <cmath>
#include <gtest/gtest.h>
#include <tuple>

using namespace magnumafcpp;

double analytic_result(double time) { return 1. / 16. * std::pow(std::pow(time, 2) + 4, 2); }

TEST(BS23, integrate_analytical_in_double) {
    const auto f = [](double t, double y0) { return t * std::sqrt(y0); };
    double t = 0;
    double y = 1;
    const double dt = 1e-2;
    const std::size_t imax = 100000;

    // initial step
    double rkerror = NAN;
    double k4 = NAN;
    std::tie(t, y, rkerror, k4) = integrator::BS23(t, dt, y, f);
    EXPECT_NEAR(y, analytic_result(t), 4e-6);

    // steps using FSAL
    for (std::size_t i = 0; i < imax; i++) {
        std::tie(t, y, rkerror, k4) = integrator::BS23(t, dt, y, k4, f);
        EXPECT_NEAR(y, analytic_result(t), 1e-0 * (i + 1) * std::pow(dt, 3)); // scaling with 3rd order in dt
        // EXPECT_NEAR(y, analytic_result(t), 4e-6); // for 1000 steps
        // std::cout << std::setprecision(32) << "i=" << i << ", t= " << t << ", y= " << y << '\n';
    }
    // Note: numerically exact same results compared to only using BS23 without k4 argument.
}

TEST(BS23, integrate_analytical_in_af_array) {
    const auto f = [](double t, const af::array& y0) { return t * af::sqrt(y0); };
    double t = 0;
    af::array y = af::constant(1., 1, f64);
    constexpr double dt = 1e-2;
    constexpr std::size_t imax = 1000;

    af::array rkerror;
    af::array k4;
    std::tie(t, y, rkerror, k4) = integrator::BS23(t, dt, y, f);
    EXPECT_NEAR(y.scalar<double>(), analytic_result(t), 4e-6);

    for (std::size_t i = 0; i < imax; i++) {
        // TODO Performance issue here:
        // call with k4 overload is far slower as without, but why? Probably memory allocation issue with af::array
        // (double version above is super fast, 31ms for 1e6 evals)
        // cpu factor 50, opencl 1.75!
        // std::tie(t, y, rkerror, k4) = integrator::BS23(t, dt, y, f); // fast (cpu: 0,117s opencl: 0,425s)
        std::tie(t, y, rkerror, k4) = integrator::BS23(t, dt, y, k4, f); // slow (cpu: 3,393s opencl: 0,742s)

        // // exactly the same numerical results as with k4 BS23:
        // const auto [t_, y_, rkerror_, k4_] = integrator::BS23(t, dt, y, k4, f);
        // t = t_;
        // y = y_;
        // rkerror = rkerror_;
        // k4 = k4_;

        EXPECT_NEAR(y.scalar<double>(), analytic_result(t), 4e-6);
        // std::cout << std::setprecision(32) << "i=" << i << ", t= " << t << ", y= " << y.scalar<double>() << '\n';
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
