#include "double_or_array.hpp"
#include <gtest/gtest.h>

using namespace magnumafcpp;

// tests commutivity of * and + operators for DOA Class, compares with DOA.get()
auto test = [](DoubleOrArray a, af::array b) {
    EXPECT_EQ((a * b).scalar<double>(), a.get(b.dims(), b.type()).scalar<double>() * b.scalar<double>());
    EXPECT_EQ((b * a).scalar<double>(), a.get(b.dims(), b.type()).scalar<double>() * b.scalar<double>());
    EXPECT_EQ((a + b).scalar<double>(), a.get(b.dims(), b.type()).scalar<double>() + b.scalar<double>());
    EXPECT_EQ((b + a).scalar<double>(), a.get(b.dims(), b.type()).scalar<double>() + b.scalar<double>());
};

TEST(util_double_or_array, double_scalar_array_op_overloads) {
    // DOA double ctor
    test(3., af::constant(2., af::dim4(1, 1, 1, 1), f64));
    test(3., af::constant(2., af::dim4(2, 2, 2, 2), f64));
    test(3., af::constant(2., af::dim4(1, 1, 1, 3), f64));

    // DOA array ctor
    test(af::constant(3., 1, 1, 1, 1, f64), af::constant(2., af::dim4(1, 1, 1, 1), f64));
    test(af::constant(3., 2, 2, 2, 2, f64), af::constant(2., af::dim4(2, 2, 2, 2), f64));
    test(af::constant(3., 1, 1, 1, 1, f64), af::constant(2., af::dim4(1, 1, 1, 3), f64));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
