#include <gtest/gtest.h>
#include "../../../src/mesh.cpp"
#include "../../../src/func.cpp"
 
// Exemplary unit test
TEST(Mesh, InitTest) {
    Mesh mesh(1,2,3,0.1,0.2,0.3);
    ASSERT_EQ(1, mesh.n0);
    ASSERT_EQ(2, mesh.n1);
    ASSERT_EQ(3, mesh.n2);
    ASSERT_EQ(0.1, mesh.dx);
    ASSERT_EQ(0.2, mesh.dy);
    ASSERT_EQ(0.3, mesh.dz);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
