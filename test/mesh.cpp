#include "mesh.hpp"
#include <gtest/gtest.h>

// Exemplary unit test
TEST(Mesh, InitTest) {
    magnumafcpp::Mesh mesh(1, 2, 3, 0.1, 0.2, 0.3);
    EXPECT_EQ(1, mesh.nx);
    EXPECT_EQ(2, mesh.ny);
    EXPECT_EQ(3, mesh.nz);
    EXPECT_EQ(0.1, mesh.dx);
    EXPECT_EQ(0.2, mesh.dy);
    EXPECT_EQ(0.3, mesh.dz);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
