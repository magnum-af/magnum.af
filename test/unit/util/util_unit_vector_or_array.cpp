#include "af_overloads.hpp"
#include "unit_vector_or_array.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace magnumafcpp;
TEST(util_vector_or_array, ctor_vec) {
    UnitVectorOrArray a(std::array<double, 3>{1, 1, 1});
    EXPECT_THAT(a.scalar_vector.value(), testing::ElementsAre(1 / std::sqrt(3), 1 / std::sqrt(3), 1 / std::sqrt(3)));
}

TEST(util_vector_or_array, ctor_array) {
    af::dim4 dims(2, 2, 2, 3);
    af::dtype type(f64);
    UnitVectorOrArray b(af::constant(1., dims, type));
    auto expected_array = af::constant(1 / std::sqrt(3), dims, type);
    EXPECT_EQ(af::sum(af::sum(af::sum(af::sum((b.array_vector.value() != expected_array).as(f64), 0), 1), 2), 3)
                  .scalar<double>(),
              0);
}

TEST(util_vector_or_array, get_as_array) {
    auto test_get_as_array = [](UnitVectorOrArray uvec, std::array<double, 3> expected_unitvec, af::dim4 dims_scalar) {
        UnitVectorOrArray v(uvec);
        auto dims_vector = dims_scalar;
        dims_vector[3] = 3;
        auto type = af::dtype(f64);
        auto testarray = v.get_as_array(dims_scalar, type);
        EXPECT_EQ(testarray.dims(), dims_vector); // Sensitive only to get_as_array for UVOR ctor with scalar vector
        auto expected_array = af::tile(af::array(1, 1, 1, 3, expected_unitvec.data()), dims_scalar);
        // check if all elements are equal:
        EXPECT_EQ(
            af::sum(af::sum(af::sum(af::sum((testarray != expected_array).as(f64), 0), 1), 2), 3).scalar<double>(), 0);
    };

    auto test_both_scalar_and_array = [test_get_as_array](std::array<double, 3> in_vec,
                                                          std::array<double, 3> expected_vec, af::dim4 dims_scalar) {
        test_get_as_array(in_vec, expected_vec, dims_scalar);
        test_get_as_array(af::tile(af::array(1, 1, 1, 3, in_vec.data()), dims_scalar), expected_vec, dims_scalar);
    };

    test_both_scalar_and_array({1, 0, 0}, {1, 0, 0}, af::dim4(1, 1, 1, 1));
    test_both_scalar_and_array({1, 0, 0}, {1, 0, 0}, af::dim4(5, 6, 7, 1));

    test_both_scalar_and_array({2, 2, 2}, {1 / std::sqrt(3), 1 / std::sqrt(3), 1 / std::sqrt(3)}, af::dim4(1, 1, 1, 1));
    test_both_scalar_and_array({2, 2, 2}, {1 / std::sqrt(3), 1 / std::sqrt(3), 1 / std::sqrt(3)}, af::dim4(2, 3, 4, 1));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
