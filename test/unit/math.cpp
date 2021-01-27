#include "math.hpp"
// #include "util/af_overloads.hpp"
#include <gtest/gtest.h>

using namespace magnumafcpp;

TEST(Math, central_diff) {
    // Checks wether central_diff of iota values (0, 1, 2, ...) is 1/h for all cells,
    // excluding boundaries.
    const auto dims = af::dim4(10, 5, 4, 6);
    const auto a = af::iota(dims).as(f64); // Note: this would somehow be f32: af::iota(dims, f64)
    const double h = 2.0;                  // step size
    const auto cent_diff = math::central_diff(a, h, math::Axis::x, math::TruncateOutput::off);
    const auto cent_diff_no_bound = cent_diff(af::seq(1, -2), af::span, af::span, af::span); // without edges
    const auto meanall = af::mean(af::mean(af::mean(af::mean(cent_diff_no_bound, 3), 2), 1), 0);
    // neglecting the boundaries, the central_diff should be 1/h at each cell
    EXPECT_EQ(meanall.scalar<double>(), 1.0 / h); // Check whether all d/dx vals are 1/h, excluding boundary
}

// Checks that af::diff1 is a forward finite difference missing the division with step size h
// i.e. diff1(f) at i is f(i+1) - f(i), not ( f(i+1) - f(i) )/h
TEST(Math, af_diff1) {
    const auto n = 10;
    auto a = af::iota(af::dim4(n, 1, 1, 1)).as(f64);
    {
        const auto diff = af::diff1(a);
        EXPECT_EQ(diff.dims(0), a.dims(0) - 1); // sesulting array has reduced dim by 1
        EXPECT_EQ(af::mean(diff).scalar<double>(), 1); // diff1 of iota is 1
    }

    // Check if af::diff1 is really a forward diff, i.e. check diff around n/2 with other value
    a(n / 2) = 1;
    {
        const auto diff_2 = af::diff1(a);
        EXPECT_EQ(diff_2(n / 2 - 1).scalar<double>(), -3.);
        EXPECT_EQ(diff_2(n / 2).scalar<double>(), 5.);
        // else it is 1:
        EXPECT_EQ(af::mean(diff_2(af::seq(0, n / 2 - 2))).scalar<double>(), 1);
        EXPECT_EQ(af::mean(diff_2(af::seq(n / 2 + 1, -1))).scalar<double>(), 1);
    }
}

TEST(Math, curl) {
    const auto dims = af::dim4(4, 5, 6, 3);
    const double d[3] = {0.1, 0.2, 0.3};
    const auto a = af::constant(1.0, dims, f64);
    const auto curl = ::magnumafcpp::math::curl_3D(a, d[0], d[1], d[2], math::TruncateOutput::off);

    // for constant field, curl should be zero (except for unhandled boundary)
    const auto curl_no_bound = curl(af::seq(1, -2), af::seq(1, -2), af::seq(1, -2), af::span); // without edges
    const auto meanall = af::mean(af::mean(af::mean(af::mean(curl_no_bound, 3), 2), 1), 0);
    EXPECT_EQ(meanall.scalar<double>(), 0.0);
    // af::print("", curl);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
