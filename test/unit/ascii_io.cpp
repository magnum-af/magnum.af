#include "ascii_io.hpp"
#include "func.hpp" // for max_abs_diff
#include <gtest/gtest.h>

using namespace magnumafcpp;

TEST(ascii_io, ascii_io_test)
{
    af::array a = af::randu(6, 5, 4, 3, f64);
    Mesh mesh(6, 5, 4, 0.1, 0.2, 0.3);

    std::string filename = ("ascii_unittest.txt");
    write_ascii(a, mesh, filename, false);

    af::array read_a;
    Mesh read_mesh(0, 0, 0, 0, 0, 0);

    read_ascii(read_a, read_mesh, filename, false);
    ASSERT_EQ(remove(filename.c_str()), 0);

    ASSERT_EQ(read_mesh.n0, 6);
    ASSERT_EQ(read_mesh.n1, 5);
    ASSERT_EQ(read_mesh.n2, 4);
    ASSERT_EQ(read_mesh.dx, 0.1);
    ASSERT_EQ(read_mesh.dy, 0.2);
    ASSERT_EQ(read_mesh.dz, 0.3);
    ASSERT_EQ(max_abs_diff(read_a, a), 0);
}


int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
