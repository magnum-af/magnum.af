#include "vtk_IO.hpp"
#include "util/func.hpp"
#include <gtest/gtest.h>
#include <stdio.h>

using namespace magnumafcpp;

TEST(vtkIO, vtiWriteReadTest) {
    af::array a = af::randu(6, 5, 4, 3, f64);
    Mesh mesh(6, 5, 4, 0.1, 0.2, 0.3);

    vti_writer_micro(a, mesh, "vti_unittest");

    af::array read_a;
    Mesh read_mesh(0, 0, 0, 0, 0, 0);

    vti_reader(read_a, read_mesh, "vti_unittest.vti");
    EXPECT_EQ(remove("vti_unittest.vti"), 0);

    EXPECT_EQ(read_mesh.nx, 6);
    EXPECT_EQ(read_mesh.ny, 5);
    EXPECT_EQ(read_mesh.nz, 4);
    EXPECT_EQ(read_mesh.dx, 0.1);
    EXPECT_EQ(read_mesh.dy, 0.2);
    EXPECT_EQ(read_mesh.dz, 0.3);
    EXPECT_EQ(max_abs_diff(read_a, a), 0);
}

TEST(vtkIO, vtrWriteReadTest) {
    af::array a = af::randu(6, 5, 4, 10, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequiMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "vtr_unittest", false);

    af::array read_a;
    NonequiMesh read_mesh(0, 0, 0, 0, {0});

    vtr_reader(read_a, read_mesh, "vtr_unittest", false);
    EXPECT_EQ(remove("vtr_unittest.vtr"), 0);

    EXPECT_EQ(read_mesh.nx, 6);
    EXPECT_EQ(read_mesh.ny, 5);
    EXPECT_EQ(read_mesh.nz, 4);
    EXPECT_EQ(read_mesh.dx, 0.1);
    EXPECT_EQ(read_mesh.dy, 0.2);

    for (unsigned i = 0; i < z_spacing.size(); i++) {
        EXPECT_NEAR(z_spacing.at(i), read_mesh.z_spacing.at(i), 2e-16);
    }

    EXPECT_EQ(max_abs_diff(read_a, a), 0);
}

TEST(vtkIO, vtrWriteReadScalarFieldTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequiMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "vtr_unittest", false);

    af::array read_a;

    vtr_reader(read_a, mesh, "vtr_unittest", false);
    EXPECT_EQ(remove("vtr_unittest.vtr"), 0);

    EXPECT_EQ(max_abs_diff(read_a, a), 0);
}

TEST(vtkIO, vtrWriteReadAddFileExtensionTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequiMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "unittesting_a", false);

    af::array read_a;

    vtr_reader(read_a, mesh, "unittesting_a.vtr", false);
    EXPECT_EQ(remove("unittesting_a.vtr"), 0);

    EXPECT_EQ(max_abs_diff(read_a, a), 0);

    vtr_writer(a, mesh, "unittesting_b.vtr", false);
    af::array read_b;
    vtr_reader(read_b, mesh, "unittesting_b", false);

    EXPECT_EQ(max_abs_diff(read_b, a), 0);
    EXPECT_EQ(remove("unittesting_b.vtr"), 0);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}