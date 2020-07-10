#include "ascii_io.hpp"
#include "func.hpp" // for max_abs_diff
#include <gtest/gtest.h>

using namespace magnumafcpp;

TEST(ascii_io, ascii_io_test_vector_field) {
    af::array a = af::randu(6, 5, 4, 3, f64);
    Mesh mesh(6, 5, 4, 0.1, 0.2, 0.3);

    std::string filename = ("ascii_unittest.txt");
    write_ascii(a, mesh, filename, false);

    af::array read_a;
    Mesh read_mesh(0, 0, 0, 0, 0, 0);

    read_ascii(read_a, read_mesh, filename, false);
    EXPECT_EQ(remove(filename.c_str()), 0);

    EXPECT_EQ(read_mesh.n0, 6);
    EXPECT_EQ(read_mesh.n1, 5);
    EXPECT_EQ(read_mesh.n2, 4);
    EXPECT_EQ(read_mesh.dx, 0.1);
    EXPECT_EQ(read_mesh.dy, 0.2);
    EXPECT_EQ(read_mesh.dz, 0.3);
    EXPECT_EQ(max_abs_diff(read_a, a), 0);
}

TEST(ascii_io, ascii_io_test_scalar_field) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    Mesh mesh(6, 5, 4, 0.1, 0.2, 0.3);

    std::string filename = ("ascii_unittest.txt");
    write_ascii(a, mesh, filename, false);

    af::array read_a;
    Mesh read_mesh(0, 0, 0, 0, 0, 0);

    read_ascii(read_a, read_mesh, filename, false);
    EXPECT_EQ(remove(filename.c_str()), 0);

    EXPECT_EQ(read_mesh.n0, 6);
    EXPECT_EQ(read_mesh.n1, 5);
    EXPECT_EQ(read_mesh.n2, 4);
    EXPECT_EQ(read_mesh.dx, 0.1);
    EXPECT_EQ(read_mesh.dy, 0.2);
    EXPECT_EQ(read_mesh.dz, 0.3);
    EXPECT_EQ(max_abs_diff(read_a, a), 0);
}

TEST(ascii_io, ascii_io_test_10d_field) {
    af::array a = af::randu(6, 5, 4, 10, f64);
    Mesh mesh(6, 5, 4, 0.2, 0.3, 0.4);

    std::string filename = ("ascii_unittest.txt");
    write_ascii(a, mesh, filename, false);

    af::array read_a;
    Mesh read_mesh(0, 0, 0, 0, 0, 0);

    read_ascii(read_a, read_mesh, filename, false);
    EXPECT_EQ(remove(filename.c_str()), 0);

    EXPECT_EQ(read_mesh.n0, 6);
    EXPECT_EQ(read_mesh.n1, 5);
    EXPECT_EQ(read_mesh.n2, 4);
    EXPECT_EQ(read_mesh.dx, 0.2);
    EXPECT_EQ(read_mesh.dy, 0.3);
    EXPECT_EQ(read_mesh.dz, 0.4);
    EXPECT_EQ(max_abs_diff(read_a, a), 0);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
