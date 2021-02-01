#include "bs23.hpp"
#include "arrayfire.h"
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
    const std::size_t imax = 1000;

    // initial step
    double rkerror;
    double k4;
    std::tie(t, y, rkerror, k4) = BS23(t, dt, y, f);
    EXPECT_NEAR(y, analytic_result(t), 4e-6);

    // steps using FSAL
    for (std::size_t i = 0; i < imax; i++) {
        std::tie(t, y, rkerror, k4) = BS23(t, dt, y, k4, f);
        EXPECT_NEAR(y, analytic_result(t), 4e-6);
        // std::cout << std::setprecision(32) << "i=" << i << ", t= " << t << ", y= " << y << '\n';
    }
    // Note: numerically exact same results compared to only using BS23 without k4 argument.
}

TEST(BS23, integrate_analytical_in_af_array) {
    const auto f = [](double t, af::array y0) { return t * af::sqrt(y0); };
    double t = 0;
    af::array y = af::constant(1., 1, f64);
    const double dt = 1e-2;
    const std::size_t imax = 1000;
    for (std::size_t i = 0; i < imax; i++) {
        af::array rkerror;
        af::array k4;
        std::tie(t, y, rkerror, k4) = BS23(t, dt, y, f);
        EXPECT_NEAR(y.scalar<double>(), analytic_result(t), 4e-6);
        // std::cout << std::setprecision(32) << "i=" << i << ", t= " << t << ", y= " << y.scalar<double>() << '\n';
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
