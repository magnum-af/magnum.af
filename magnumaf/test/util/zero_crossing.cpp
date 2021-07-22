#include "util/zero_crossing.hpp"
#include <gtest/gtest.h>

using namespace magnumaf;

double x(double x) { return x; }

TEST(ZeroCrossing, x) {
    double precision = 1e-12;
    ZeroCrossing zc(x, precision, 10, -1, 1.01, 100, false);
    auto result = zc.calc_x_and_f();
    EXPECT_NEAR(0, result.first, precision);
    EXPECT_NEAR(0, result.second, precision);
}

TEST(ZeroCrossing, x_symmetric) {
    double precision = 1e-12;
    ZeroCrossing zc(x, precision, 10, -1, 1, 100, false);
    auto result = zc.calc_x_and_f();
    EXPECT_NEAR(0, result.first, precision);
    EXPECT_NEAR(0, result.second, precision);
}

double x_plus_1(double x) { return x + 1; }

TEST(ZeroCrossing, x_plus_1) {
    double precision = 1e-12;
    ZeroCrossing zc(x_plus_1, precision, 10, -2.1, 2, 100, false);
    auto result = zc.calc_x_and_f();
    EXPECT_NEAR(-1, result.first, precision);
    EXPECT_NEAR(0, result.second, precision);
}

TEST(ZeroCrossing, x_plus_1_negative_range) {
    double precision = 1e-12;
    ZeroCrossing zc(x_plus_1, precision, 10, 2.1, -2, 100, false);
    auto result = zc.calc_x_and_f();
    EXPECT_NEAR(-1, result.first, precision);
    EXPECT_NEAR(0, result.second, precision);
}

double minus_x_plus_1(double x) { return -(x + 1); }

TEST(ZeroCrossing, minus_x_plus_1) {
    double precision = 1e-12;
    ZeroCrossing zc(minus_x_plus_1, precision, 10, -2.1, 2, 100, false);
    auto result = zc.calc_x_and_f();
    EXPECT_NEAR(-1, result.first, precision);
    EXPECT_NEAR(0, result.second, precision);
}

// f(x) = 0, get from range [1,2] to range covering 0
TEST(ZeroCrossing, x_out_of_range) {
    double precision = 1e-12;
    ZeroCrossing zc(x, precision, 10, 1, 2, 100, false);
    auto result = zc.calc_x_and_f();
    EXPECT_NEAR(0, result.first, precision);
    EXPECT_NEAR(0, result.second, precision);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
