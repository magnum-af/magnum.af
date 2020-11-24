#include "field_terms/micro/demag_field.hpp"
#include <gtest/gtest.h>

using namespace magnumafcpp;

TEST(MicroDemag, EnergyOfHomogeneousCube) {
    const double x = 1.e-9, y = 1.e-9, z = 1.e-9;
    const int nx = 8, ny = 10, nz = 12;

    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
    double Ms = 1e5;

    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(af::span, af::span, af::span, 0) = 1;
    State state(mesh, Ms, m);
    DemagField demag(mesh);

    std::cout.precision(24);
    double demagE = demag.E(state);
    double analytic = 1. / 6. * x * y * z * pow(state.Ms, 2) * constants::mu0;

    // EXPECT_NEAR(llgE, analytic, (llgE+analytic)/2. * 1e-12);// cpu precision
    EXPECT_NEAR(demagE, analytic,
                (demagE + analytic) / 2. * 5.3e-8); // opencl precision

    double demagEh = demag.E(state, demag.h(state));
    EXPECT_NEAR(demagEh, analytic,
                (demagEh + analytic) / 2. * 5.3e-8); // opencl precision
}

TEST(MicroDemag, EnergyOfHomogeneousCubeWithAirbox) {
    const double x = 1.e-9, y = 1.e-9, z = 1.e-9;
    const int nx = 12, ny = 10, nz = 8;

    Mesh mesh(2 * nx, 2 * ny, 2 * nz, x / nx, y / ny, z / nz);
    double Ms = 1e5;

    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(af::seq(nx / 2, 2 * nx - nx / 2 - 1), af::seq(ny / 2, 2 * ny - ny / 2 - 1), af::seq(nz / 2, 2 * nz - nz / 2 - 1),
      0) = 1;
    State state(mesh, Ms, m, false);
    DemagField demag(mesh);

    std::cout.precision(24);
    double demagE = demag.E(state);
    double analytic = 1. / 6. * x * y * z * pow(state.Ms, 2) * constants::mu0;

    // EXPECT_NEAR(llgE, analytic, (llgE+analytic)/2. * 1e-12);// cpu precision
    EXPECT_NEAR(demagE, analytic,
                (demagE + analytic) / 2. * 5.3e-8); // opencl precision

    double demagEh = demag.E(state, demag.h(state));
    EXPECT_NEAR(demagEh, analytic,
                (demagEh + analytic) / 2. * 5.3e-8); // opencl precision
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
