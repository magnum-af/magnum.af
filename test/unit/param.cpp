#include <gtest/gtest.h>
#include <math.h>
#include "../../src/material.hpp"
#include "../../src/constants.hpp"

using namespace magnumafcpp;

TEST(Material, Init_p_Test)
{
    Material material;
    material.p = 2.5;
    ASSERT_EQ(material.p, 2.5);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
