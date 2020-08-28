#include "cubic_anisotropy_field.hpp"
#include "util.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
using namespace magnumafcpp;
TEST(CubicAnisotropyField, Constructor) {
    CubicAnisotropyField caniso(1, 1, 1, {1, 1, 0}, {-1, 1, 0});
    EXPECT_THAT(caniso.c1, testing::ElementsAre(1 / std::sqrt(2), 1 / std::sqrt(2), 0));
    EXPECT_THAT(caniso.c2, testing::ElementsAre(-1 / std::sqrt(2), 1 / std::sqrt(2), 0));
    EXPECT_THAT(caniso.c3, testing::ElementsAre(0, 0, 1));
}

void Kc1_test(const double Kc1) {
    const double Kc2 = 0;
    const double Kc3 = 0;
    CubicAnisotropyField caniso(Kc1, Kc2, Kc3, {1, 0, 0}, {0, 1, 0});
    const Mesh mesh(1, 1, 1, 1e-9, 1e-9, 1e-9);
    const double Ms = 1 / constants::mu0;
    af::array m = af::constant(0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(mesh, Ms, m);

    // Spin m in xy-plane
    const unsigned imax = 64;
    for (unsigned i = 0; i < imax; ++i) {
        double x = std::sin(2 * M_PI * i / imax);
        double y = std::cos(2 * M_PI * i / imax);
        state.m(0, 0, 0, 0) = x;
        state.m(0, 0, 0, 1) = y;

        std::array<double, 3> m = {x, y, 0};
        double c1m2 = std::pow(dot_product(caniso.c1, m), 2);
        double c2m2 = std::pow(dot_product(caniso.c2, m), 2);
        double c3m2 = std::pow(dot_product(caniso.c3, m), 2);

        double E_density_analytic = caniso.Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2);
        double E_analytic = E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
        double E_calculated = caniso.E(state);
        double E_density_calculated = E_calculated / (mesh.dx * mesh.dy * mesh.dz);
        EXPECT_NEAR(E_calculated, E_analytic, 1e-37);
        EXPECT_NEAR(E_density_calculated, E_density_analytic, 1e-10);
        // std::cout << i << " " << E_density_analytic << " " << E_density_calculated << std::endl;
        // std::cout << i << " " << E_analytic << " " << E_calculated << std::endl << std::endl;
    }
}

TEST(CubicAnisotropyField, Kc1_energy_test) {
    Kc1_test(1e6);
    Kc1_test(-1e6);
}

void Kc3_test(const double Kc3) {
    const double Kc1 = 0;
    const double Kc2 = 0;
    const Mesh mesh(1, 1, 1, 1e-9, 1e-9, 1e-9);
    const double Ms = 1 / constants::mu0;
    af::array m = af::constant(0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(mesh, Ms, m);

    // Spin m in xy-plane
    const unsigned imax = 64;
    for (unsigned i = 0; i < imax; ++i) {
        double x = std::sin(2 * M_PI * i / imax);
        double y = std::cos(2 * M_PI * i / imax);
        state.m(0, 0, 0, 0) = x;
        state.m(0, 0, 0, 1) = y;

        std::array<double, 3> m = {x, y, 0};
        CubicAnisotropyField caniso(Kc1, Kc2, Kc3, {1, 0, 0}, {0, 1, 0});
        double c1m4 = std::pow(dot_product(caniso.c1, m), 4);
        double c2m4 = std::pow(dot_product(caniso.c2, m), 4);
        double c3m4 = std::pow(dot_product(caniso.c3, m), 4);

        double Kc3_E_density_analytic = caniso.Kc3 * (c1m4 * c2m4 + c1m4 * c3m4 + c2m4 * c3m4);
        double Kc3_E_analytic = Kc3_E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
        double Kc3_E_calculated = caniso.E(state);
        double Kc3_E_density_calculated = Kc3_E_calculated / (mesh.dx * mesh.dy * mesh.dz);
        EXPECT_NEAR(Kc3_E_calculated, Kc3_E_analytic, 1e-37);
        EXPECT_NEAR(Kc3_E_density_calculated, Kc3_E_density_analytic, 1e-10);
        // std::cout << i << " " << Kc3_E_density_analytic << " " << Kc3_E_density_calculated << std::endl;
        // std::cout << i << " " << Kc3_E_analytic << " " << Kc3_E_calculated << std::endl << std::endl;
    }
}

TEST(CubicAnisotropyField, Kc3_energy_test) {
    Kc3_test(1e6);
    Kc3_test(-1e6);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
