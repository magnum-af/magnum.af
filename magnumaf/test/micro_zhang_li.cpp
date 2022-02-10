#include "field_terms/micro/spin_transfer_torque_zhang_li_field.hpp"
#include "mesh.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace magnumaf;

TEST(MicroZhangLi, HeffTest) {
    const double x = 1.e-9, y = 1.e-9, z = 1.e-9;
    // const int nx = 8, ny = 10, nz = 12;
    const int nx = 4, ny = 3, nz = 2;

    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
    // double Ms = 1e5;
    // double Ms = 0;
    af::array Ms = af::constant(0, mesh::dims_s(mesh), f64);
    Ms(af::seq(0, nx / 2), af::span, af::span) = 1e5;

    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(af::span, af::span, af::span, 0) = 1;
    State state(mesh, Ms, m);

    af::array j = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    j(af::span, af::span, af::span, 0) = 1;
    double beta = 1.0;
    double xi = 1.0;
    auto zhang_li = SpinTransferTorqueZhangLiField(j, beta, xi);
    auto h = zhang_li.H_in_Apm(state);
    print("h", h);
    std::vector<double> h_expected = {2., 2., 2.};
    EXPECT_NEAR(h(0, 0, 0, 0).scalar<double>(), h_expected[0], 1e-5);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
