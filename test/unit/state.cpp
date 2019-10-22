#include "../../src/mesh.cpp"
#include "../../src/nonequispaced_mesh.cpp"
#include "../../src/func.cpp"
#include "../../src/state.cpp"
#include "../../src/vtk_IO.cpp"
//#include "arrayfire.h"
#include <gtest/gtest.h>
#include <vector>

using namespace magnumafcpp;

TEST(State, integral_nonequimesh_dz_as_gauss_sum) {
    const double dx=1e-9, dy=2e-9, dz=3e-9;
    const int nx=1, ny=1, nz=10;
    std::vector<double> z_spacing;
    for(int i=1; i<=nz; i++){
        z_spacing.push_back(i * dz);// gauss sum: dz * (1 + 2 + 3 + ... + n)
    }
    NonequispacedMesh mesh(nx, ny, dx, dy, z_spacing);
    af::array m = af::constant(0.0, mesh.dims, f64);
    for(int i=0; i<nz; i++){
        if( i%2 == 0){
            m(af::span, af::span, i, 2) = 1;
        }
        else {
            m(af::span, af::span, i, 2) = -1;
        }
    }
    State state(mesh, 1., m);
    //calculating gauss sums s_n = sum_{i=1}^n (a_1 + (i-1)d) = n \frac{a_1 + a_n}{2}
    double positive_sum =   nz/2 * (1 + (nz-1))/2;// 1 + 3 + 5 ... + 9
    double negative_sum = - nz/2 * (2 + nz)/2;// 2 + 4 + 6 ... + 10
    double result = dx * nx * dy * ny * dz * (positive_sum + negative_sum);
    ASSERT_NEAR(state.integral_nonequimesh(m), result, 1e-35);
}


TEST(State, vtkIO_vtrWriteReadTest) {
    af::array a = af::randu(6, 5, 4, 10, f64);
    std::vector<double> z_spacing = {0.1, 0.2, 0.3, 0.4};
    NonequispacedMesh mesh(6, 5, 0.1, 0.2, z_spacing);
    State state(mesh, 0, a, false, true);

    state.vtr_writer("vtr_unittest");

    State read_state(NonequispacedMesh(0, 0, 0, 0, {0}), 4e5, af::constant(0, 1, f64), false, true);
    read_state.vtr_reader("vtr_unittest");
    EXPECT_EQ(remove( "vtr_unittest.vtr" ), 0);

    EXPECT_EQ(read_state.nonequimesh.nx, 6);
    EXPECT_EQ(read_state.nonequimesh.ny, 5);
    EXPECT_EQ(read_state.nonequimesh.nz, 4);
    EXPECT_EQ(read_state.nonequimesh.dx, 0.1);
    EXPECT_EQ(read_state.nonequimesh.dy, 0.2);

    for (unsigned i = 0; i < z_spacing.size(); i++){
        ASSERT_NEAR(z_spacing.at(i), read_state.nonequimesh.z_spacing.at(i), 2e-16);
    }

    ASSERT_EQ(max_abs_diff(read_state.m, a), 0);
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
