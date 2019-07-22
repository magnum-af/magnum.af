#include "../../../src/mesh.cpp"
#include "../../../src/nonequispaced_mesh.cpp"
#include "../../../src/func.cpp"
#include "../../../src/state.cpp"
#include "../../../src/vtk_IO.cpp"
//#include "arrayfire.h"
#include <gtest/gtest.h>
#include <vector>

using namespace magnumaf;

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

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
