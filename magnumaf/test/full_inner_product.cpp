#include "arrayfire.h"
#include "util/util.hpp"
#include <gtest/gtest.h>

using namespace magnumaf;

// Exemplary unit test
TEST(FullInnerProductTest, n) {
    af::array a = af::constant(1., 3, 3, 3, 3, f64);
    af::array b = af::constant(2., 3, 3, 3, 3, f64);
    EXPECT_EQ(3 * 3 * 3 * 3 * 2., util::full_inner_product(a, b));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
