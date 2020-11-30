#include "double_or_array.hpp"
#include <gtest/gtest.h>

using namespace magnumafcpp;

// tests commutivity of * and + operators for DOA Class, compares with DOA.get()
auto test = [](DoubleOrArray a, af::array b) {
    const double aa = a.get(b.dims(), b.type()).scalar<double>();
    const double bb = b.scalar<double>();
    EXPECT_EQ((a + b).scalar<double>(), aa + bb);
    EXPECT_EQ((b + a).scalar<double>(), bb + aa);
    EXPECT_EQ((a - b).scalar<double>(), aa - bb);
    // Note: test if enabled:// EXPECT_EQ((b - a).scalar<double>(), bb - aa);
    EXPECT_EQ((a * b).scalar<double>(), aa * bb);
    EXPECT_EQ((b * a).scalar<double>(), bb * aa);
    EXPECT_EQ((a / b).scalar<double>(), aa / bb);
    // Note: test if enabled:// EXPECT_EQ((b / a).scalar<double>(), bb / aa);
};

TEST(util_double_or_array, double_scalar_array_op_overloads) {
    // DOA double ctor
    test(DoubleOrArray{3.}, af::constant(2., af::dim4(1, 1, 1, 1), f64));
    test(DoubleOrArray{3.}, af::constant(2., af::dim4(1, 1, 1, 3), f64));
    // this works when ctor is not explicit with implicit convert via ctor
    // test(3., af::constant(2., af::dim4(1, 1, 1, 1), f64));
    // test(3., af::constant(2., af::dim4(1, 1, 1, 3), f64));

    // DOA array ctor
    test(DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(1, 1, 1, 1), f64));
    test(DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(1, 1, 1, 3), f64));
    // this works when ctor is not explicit with implicit convert via ctor
    // test(af::constant(3., 1, 1, 1, 1, f64), af::constant(2., af::dim4(1, 1, 1, 1), f64));
    // test(af::constant(3., 1, 1, 1, 1, f64), af::constant(2., af::dim4(1, 1, 1, 3), f64));

    //// checking exceptions in ctor
    // test(DoubleOrArray{3.}, af::constant(2., af::dim4(2, 2, 2, 2), f64));
    // test(DoubleOrArray{af::constant(3., 2, 2, 2, 2, f64)}, af::constant(2., af::dim4(2, 2, 2, 2), f64));
    //// as well runtime
    // test(DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(2, 2, 2, 1), f64));
    // test(DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(1, 1, 1, 4), f64));
    //// gtest ustream this could be used:
    // EXPECT_THAT([]() { throw std::runtime_error("message"); }, Throws<std::runtime_error>());
}

// Demo for unwanted implicit conversion from double to af::array
// via conversion ctor af::array(dim_t dim0, dtype ty=f32)
// Resolved by deleting all operators for operator(double, DoubleOrArray)
// Should be resolved by declaring af::array() ctor explicit

// TEST(util_double_or_array, error_double_to_array_promotion_demo) {
//    const double d = 5.;
//    // delete op affects int also:// const int d = 5.;
//    const DoubleOrArray a{2.};
//    af::print("", d + a);
//    af::print("", d * a);
//    // af::print("", d - a);
//    // af::print("", d / a);
//
//    af::print("", a + d);
//    af::print("", a * d);
//    af::print("", a - d);
//    af::print("", a / d);
//}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
