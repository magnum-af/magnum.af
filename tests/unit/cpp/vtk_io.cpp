#include <gtest/gtest.h>
#include "../../../src/mesh.cpp"
#include "../../../src/func.cpp"
#include "../../../src/vtk_IO.cpp"
 
TEST(vtkIO, vtrWriteReadTest) {
    af::array a = af::randu(6, 5, 4, 10, f64);
    Mesh mesh(6, 5, 4, 0.1, 0.2, 0);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};

    vtr_writer(a, mesh, z_spacing, "vtr_unittest", false);

    af::array read_a;
    Mesh read_mesh(0, 0, 0, 0, 0, 0);
    std::vector<double> read_z_spacing;

    vtr_reader(read_a, read_mesh, read_z_spacing, "vtr_unittest", false);

    ASSERT_EQ(read_mesh.n0, 6);
    ASSERT_EQ(read_mesh.n1, 5);
    ASSERT_EQ(read_mesh.n2, 4);
    ASSERT_EQ(read_mesh.dx, 0.1);
    ASSERT_EQ(read_mesh.dy, 0.2);
    ASSERT_EQ(read_mesh.dz, 0);

    for (unsigned i = 0; i < z_spacing.size(); i++){
        ASSERT_NEAR(z_spacing.at(i), read_z_spacing.at(i), 2e-16); 
    }

    ASSERT_EQ(max_abs_diff(read_a, a), 0);
}


TEST(vtkIO, vtrWriteReadScalarFieldTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    Mesh mesh(6, 5, 4, 0.1, 0.2, 0);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};

    vtr_writer(a, mesh, z_spacing, "vtr_unittest", false);

    af::array read_a;

    vtr_reader(read_a, mesh, z_spacing, "vtr_unittest", false);

    ASSERT_EQ(max_abs_diff(read_a, a), 0);
}


TEST(vtkIO, vtrWriteReadAddFileExtensionTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    Mesh mesh(6, 5, 4, 0.1, 0.2, 0);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};

    vtr_writer(a, mesh, z_spacing, "unittesting_a", false);

    af::array read_a;

    vtr_reader(read_a, mesh, z_spacing, "unittesting_a.vtr", false);

    ASSERT_EQ(max_abs_diff(read_a, a), 0);

    vtr_writer(a, mesh, z_spacing, "unittesting_b.vtr", false);
    af::array read_b;
    vtr_reader(read_b, mesh, z_spacing, "unittesting_b", false);

    ASSERT_EQ(max_abs_diff(read_b, a), 0);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
