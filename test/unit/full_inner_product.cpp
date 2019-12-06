#include <gtest/gtest.h>
#include "arrayfire.h"
#include "../../src/func.hpp"

using namespace magnumafcpp;


// Exemplary unit test
TEST(FullInnerProductTest, n) {
    af::array a = af::constant(1., 3, 3, 3, 3, f64);
    af::array b = af::constant(2., 3, 3, 3, 3, f64);
    ASSERT_EQ(3*3*3*3*2., full_inner_product(a, b));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
