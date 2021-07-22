#include "rkf45.hpp"
#include "arrayfire.h"
#include "math.h"

#include <cmath>
#include <gtest/gtest.h>
#include <tuple>

using namespace magnumaf;

double analytic_result(double time) { return 1. / 16. * std::pow(std::pow(time, 2) + 4, 2); }

TEST(RKF45, integrate_analytical_in_double) {
    auto f = [](double t, double y0) { return t * std::sqrt(y0); };
    double t = 0;
    double y = 1;
    const double dt = 1e-2;
    const std::size_t imax = 1000;
    for (std::size_t i = 0; i < imax; i++) {
        double rkerror = NAN;
        std::tie(t, y, rkerror) = RKF45(t, dt, y, f);
        EXPECT_NEAR(y, analytic_result(t), 1e-8);
        // std::cout << std::setprecision(32) << "i=" << i << ", t= " << t << ", y= " << y << '\n';
    }
}

TEST(RKF45, integrate_analytical_in_af_array) {
    auto f = [](double t, const af::array& y0) { return t * af::sqrt(y0); };
    double t = 0;
    af::array y = af::constant(1., 1, f64);
    const double dt = 1e-2;
    const std::size_t imax = 100;
    for (std::size_t i = 0; i < imax; i++) {
        af::array rkerror;
        std::tie(t, y, rkerror) = RKF45(t, dt, y, f);
        EXPECT_NEAR(y.scalar<double>(), analytic_result(t), 1e-8);
        // std::cout << std::setprecision(32) << "i=" << i << ", t= " << t << ", y= " << y.scalar<double>() << '\n';
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
