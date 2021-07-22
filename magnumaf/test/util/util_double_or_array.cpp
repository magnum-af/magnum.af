#include "double_or_array.hpp"
#include <gtest/gtest.h>

using namespace magnumaf;

// tests commutivity of * and + operators for DOA Class, compares with DOA()
auto test = [](const util::DoubleOrArray& a, const af::array& b) {
    const auto aa = a(b.dims(), b.type()).scalar<double>();
    const auto bb = b.scalar<double>();
    EXPECT_EQ((a + b).scalar<double>(), aa + bb);
    EXPECT_EQ((b + a).scalar<double>(), bb + aa);
    EXPECT_EQ((a - b).scalar<double>(), aa - bb);
    EXPECT_EQ((b - a).scalar<double>(), bb - aa);
    EXPECT_EQ((a * b).scalar<double>(), aa * bb);
    EXPECT_EQ((b * a).scalar<double>(), bb * aa);
    EXPECT_EQ((a / b).scalar<double>(), aa / bb);
    EXPECT_EQ((b / a).scalar<double>(), bb / aa);
};

TEST(util_double_or_array, double_scalar_array_op_overloads) {
    // DOA double ctor
    test(util::DoubleOrArray{3.}, af::constant(2., af::dim4(1, 1, 1, 1), f64));
    test(util::DoubleOrArray{3.}, af::constant(2., af::dim4(1, 1, 1, 3), f64));
    // this works when ctor is not explicit with implicit convert via ctor
    // test(3., af::constant(2., af::dim4(1, 1, 1, 1), f64));
    // test(3., af::constant(2., af::dim4(1, 1, 1, 3), f64));

    // DOA array ctor
    test(util::DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(1, 1, 1, 1), f64));
    test(util::DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(1, 1, 1, 3), f64));
    // this works when ctor is not explicit with implicit convert via ctor
    // test(af::constant(3., 1, 1, 1, 1, f64), af::constant(2., af::dim4(1, 1, 1, 1), f64));
    // test(af::constant(3., 1, 1, 1, 1, f64), af::constant(2., af::dim4(1, 1, 1, 3), f64));

    // testing af::array operator/(const af::array& a, const DoubleOrArray& b):
    {
        const int nx = 5, ny = 4, nz = 3;
        const auto dims = af::dim4(nx, ny, nz, 1);
        af::array Ms_init = af::constant(2., dims, f64);
        Ms_init(1, af::span, af::span) = 0.0;
        const auto Ms = util::DoubleOrArray{Ms_init};
        const auto a = af::constant(2, dims, f64);
        test(Ms, a);
        const auto allsum = [](const af::array& a) { return af::sum(af::sum(af::sum(af::sum(a, 0), 1), 2), 3); };
        // af::print("test", a / Ms);
        EXPECT_EQ(allsum(a / Ms).scalar<double>(), (nx - 1) * ny * nz);
    }

    //// checking exceptions in ctor
    // test(util::DoubleOrArray{3.}, af::constant(2., af::dim4(2, 2, 2, 2), f64));
    // test(util::DoubleOrArray{af::constant(3., 2, 2, 2, 2, f64)}, af::constant(2., af::dim4(2, 2, 2, 2), f64));
    //// as well runtime
    // test(util::DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(2, 2, 2, 1), f64));
    // test(util::DoubleOrArray{af::constant(3., 1, 1, 1, 1, f64)}, af::constant(2., af::dim4(1, 1, 1, 4), f64));
    //// gtest ustream this could be used:
    // EXPECT_THAT([]() { throw std::runtime_error("message"); }, Throws<std::runtime_error>());
}

// Demo for unwanted implicit conversion from double to af::array
// via conversion ctor af::array(dim_t dim0, dtype ty=f32)
// Resolved by deleting all operators for operator(double, util::DoubleOrArray)
// Should be resolved by declaring af::array() ctor explicit

// TEST(util_double_or_array, error_double_to_array_promotion_demo) {
//    const double d = 5.;
//    // delete op affects int also:// const int d = 5.;
//    const util::DoubleOrArray a{2.};
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
