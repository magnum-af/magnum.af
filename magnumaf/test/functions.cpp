#include "math.hpp"
#include "util/util.hpp" // legacy location partial_mean
#include <array>
#include <gtest/gtest.h>

using namespace magnumaf;

TEST(Func, Cross4) {
    const int x = 1, y = 1, z = 1;
    const std::array<double, 3> aval = {1, 2, 3};
    const std::array<double, 3> bval = {4, 5, 6};
    af::array a(x, y, z, 3, aval.data());
    af::array b(x, y, z, 3, bval.data());
    af::array c = math::cross4(a, b);
    EXPECT_EQ(c(0, 0, 0, 0).scalar<double>(), -3.);
    EXPECT_EQ(c(0, 0, 0, 1).scalar<double>(), 6.);
    EXPECT_EQ(c(0, 0, 0, 2).scalar<double>(), -3.);
}

TEST(Func, Cross4shift) {
    const int x = 1, y = 1, z = 1;
    const std::array<double, 3> aval = {1, 2, 3};
    const std::array<double, 3> bval = {4, 5, 6};
    af::array a(x, y, z, 3, aval.data());
    af::array b(x, y, z, 3, bval.data());
    af::array c = math::cross4shift(a, b);
    EXPECT_EQ(c(0, 0, 0, 0).scalar<double>(), -3.);
    EXPECT_EQ(c(0, 0, 0, 1).scalar<double>(), 6.);
    EXPECT_EQ(c(0, 0, 0, 2).scalar<double>(), -3.);
}

TEST(Func, Cross4shift_nxzy) {
    const std::array<double, 3> aval = {1, 2, 3};
    const std::array<double, 3> bval = {4, 5, 6};
    af::array a(1, 1, 1, 3, aval.data());
    af::array b(1, 1, 1, 3, bval.data());
    const int x = 4, y = 8, z = 9;
    a = af::tile(a, x, y, z, 1);
    b = af::tile(b, x, y, z, 1);
    af::array c = math::cross4shift(a, b);
    EXPECT_EQ(af::mean(af::mean(af::mean(c(0, 0, 0, 0), 0), 1), 2).scalar<double>(), -3.);
    EXPECT_EQ(af::mean(af::mean(af::mean(c(0, 0, 0, 1), 0), 1), 2).scalar<double>(), 6.);
    EXPECT_EQ(af::mean(af::mean(af::mean(c(0, 0, 0, 2), 0), 1), 2).scalar<double>(), -3.);
}

TEST(Func, spacial_mean_in_region) {
    const auto test = [](af::dtype vec_type, af::dtype region_type = af::dtype::u32) {
        af::array vectorfield = af::constant(0.0, 6, 1, 1, 3, vec_type);
        vectorfield(0, 0, 0, 0) = 2.5;
        vectorfield(1, 0, 0, 0) = 2.5;
        vectorfield(2, 0, 0, 0) = 2.5;

        vectorfield(0, 0, 0, 1) = 3.5;
        vectorfield(1, 0, 0, 1) = 3.5;
        vectorfield(2, 0, 0, 1) = 3.5;

        vectorfield(0, 0, 0, 2) = 4.5;
        vectorfield(1, 0, 0, 2) = 4.5;
        vectorfield(2, 0, 0, 2) = 4.5;

        af::array region = af::constant(0, 6, 1, 1, 1, region_type);
        region(0) = 1;
        region(1) = 2; // only non-zero is checked: non-one values are converted to
                       // ones, so this is also fine
        region(2) = 3;
        auto mean = magnumaf::util::spacial_mean_in_region(vectorfield, region);
        EXPECT_EQ(mean[0], 2.5);
        EXPECT_EQ(mean[1], 3.5);
        EXPECT_EQ(mean[2], 4.5);
    };
    test(af::dtype::f64);
    test(af::dtype::f32);

    if (af::isHalfAvailable(af::getDevice())) {
        test(af::dtype::f16);
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
