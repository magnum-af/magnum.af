#include "vtk_io.hpp"
#include "util/util.hpp"
#include <chrono> // for sleep
#include <cstdio>
#include <gtest/gtest.h>
#include <thread> // for sleep
// #include <future> // for async

using namespace magnumaf;

// TODO sleep is bad, enable with std::async call using future to wait for filewrite
// TEST(vtkIO, async_vtiWriteReadTest) {
//     const auto dims = af::dim4(60, 5, 4, 3);
//     const auto sleep_milis = 500; // for  60, 5, 4, 3, async takes 0.005614s
//     const double dx = 0.1, dy = 0.2, dz = 0.3;
//     af::array a = af::randu(dims, f64);
//     const auto a_ref = a;
//     Mesh mesh(dims.dims[0], dims.dims[1], dims.dims[2], dx, dy, dz);
//
//     // auto timer = af::timer::start();
//     async_vti_writer_micro(a, mesh, "vti_unittest.vti");
//     // std::cout << "async vti [s]: " << timer.stop() << std::endl;
//
//     a = 0; // Note: modifying a while async write runs must not change the output
//     std::this_thread::sleep_for(std::chrono::milliseconds(sleep_milis)); // waiting for output to finish
//
//     // // Alternative: use future and wait for it:
//     // auto future = std::async(std::launch::async, vti_writer_micro, a, mesh, "vti_unittest.vti");
//     // future.wait();
//
//     af::array read_a;
//     Mesh read_mesh(0, 0, 0, 0, 0, 0);
//
//     vti_reader(read_a, read_mesh, "vti_unittest.vti");
//     EXPECT_EQ(remove("vti_unittest.vti"), 0);
//
//     EXPECT_EQ(read_mesh.nx, dims.dims[0]);
//     EXPECT_EQ(read_mesh.ny, dims.dims[1]);
//     EXPECT_EQ(read_mesh.nz, dims.dims[2]);
//     EXPECT_EQ(read_mesh.dx, dx);
//     EXPECT_EQ(read_mesh.dy, dy);
//     EXPECT_EQ(read_mesh.dz, dz);
//     EXPECT_EQ(util::max_abs_diff(read_a, a_ref), 0);
// }

TEST(vtkIO, vtiWriteReadTest) {
    const auto dims = af::dim4(6, 5, 4, 3);
    const double dx = 0.1, dy = 0.2, dz = 0.3;
    af::array a = af::randu(dims, f64);
    Mesh mesh(dims.dims[0], dims.dims[1], dims.dims[2], dx, dy, dz);

    vti_writer_micro(a, mesh, "vti_unittest");

    af::array read_a;
    Mesh read_mesh(0, 0, 0, 0, 0, 0);

    vti_reader(read_a, read_mesh, "vti_unittest.vti");
    EXPECT_EQ(remove("vti_unittest.vti"), 0);

    EXPECT_EQ(read_mesh.nx, dims.dims[0]);
    EXPECT_EQ(read_mesh.ny, dims.dims[1]);
    EXPECT_EQ(read_mesh.nz, dims.dims[2]);
    EXPECT_EQ(read_mesh.dx, dx);
    EXPECT_EQ(read_mesh.dy, dy);
    EXPECT_EQ(read_mesh.dz, dz);
    EXPECT_EQ(util::max_abs_diff(read_a, a), 0);
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

    EXPECT_EQ(util::max_abs_diff(read_a, a), 0);
}

TEST(vtkIO, vtrWriteReadScalarFieldTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequiMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "vtr_unittest", false);

    af::array read_a;

    vtr_reader(read_a, mesh, "vtr_unittest", false);
    EXPECT_EQ(remove("vtr_unittest.vtr"), 0);

    EXPECT_EQ(util::max_abs_diff(read_a, a), 0);
}

TEST(vtkIO, vtrWriteReadAddFileExtensionTest) {
    af::array a = af::randu(6, 5, 4, 1, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequiMesh mesh(6, 5, 0.1, 0.2, z_spacing);

    vtr_writer(a, mesh, "unittesting_a", false);

    af::array read_a;

    vtr_reader(read_a, mesh, "unittesting_a.vtr", false);
    EXPECT_EQ(remove("unittesting_a.vtr"), 0);

    EXPECT_EQ(util::max_abs_diff(read_a, a), 0);

    vtr_writer(a, mesh, "unittesting_b.vtr", false);
    af::array read_b;
    vtr_reader(read_b, mesh, "unittesting_b", false);

    EXPECT_EQ(util::max_abs_diff(read_b, a), 0);
    EXPECT_EQ(remove("unittesting_b.vtr"), 0);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
