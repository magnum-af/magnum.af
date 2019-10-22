#include "../../../src/mesh.cpp"
#include "../../../src/nonequispaced_mesh.cpp"
#include "../../../src/func.cpp"
#include "../../../src/vtk_IO.cpp"
#include <gtest/gtest.h>
#include <stdio.h>

using namespace magnumafcpp;


TEST(vtkIO, vtrWriteReadTest) {
    af::array a = af::randu(6, 5, 4, 10, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequispacedMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "vtr_unittest", false);

    af::array read_a;
    NonequispacedMesh read_mesh(0, 0, 0, 0, {0});

    vtr_reader(read_a, read_mesh, "vtr_unittest", false);
    ASSERT_EQ(remove( "vtr_unittest.vtr" ), 0);

    ASSERT_EQ(read_mesh.nx, 6);
    ASSERT_EQ(read_mesh.ny, 5);
    ASSERT_EQ(read_mesh.nz, 4);
    ASSERT_EQ(read_mesh.dx, 0.1);
    ASSERT_EQ(read_mesh.dy, 0.2);

    for (unsigned i = 0; i < z_spacing.size(); i++){
        ASSERT_NEAR(z_spacing.at(i), read_mesh.z_spacing.at(i), 2e-16);
    }

    ASSERT_EQ(max_abs_diff(read_a, a), 0);
}


TEST(vtkIO, vtrWriteReadScalarFieldTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequispacedMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "vtr_unittest", false);

    af::array read_a;

    vtr_reader(read_a, mesh, "vtr_unittest", false);
    ASSERT_EQ(remove( "vtr_unittest.vtr" ), 0);

    ASSERT_EQ(max_abs_diff(read_a, a), 0);
}


TEST(vtkIO, vtrWriteReadAddFileExtensionTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequispacedMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "unittesting_a", false);

    af::array read_a;

    vtr_reader(read_a, mesh, "unittesting_a.vtr", false);
    ASSERT_EQ(remove( "unittesting_a.vtr" ), 0);

    ASSERT_EQ(max_abs_diff(read_a, a), 0);

    vtr_writer(a, mesh, "unittesting_b.vtr", false);
    af::array read_b;
    vtr_reader(read_b, mesh, "unittesting_b", false);

    ASSERT_EQ(max_abs_diff(read_b, a), 0);
    ASSERT_EQ(remove( "unittesting_b.vtr" ), 0);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
