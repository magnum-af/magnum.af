#include <gtest/gtest.h>
#include <math.h>
#include "../../../src/param.cpp"
 
// Exemplary unit test
TEST(ParamInitTest, mu0) {
    Param param;
    ASSERT_EQ(4e-7 * M_PI, param.mu0);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
