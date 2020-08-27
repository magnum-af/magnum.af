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

// double E_analytic_Ku1(unsigned i, double Ku1){
//        double x = std::sin(2 * M_PI * i);
//        double y = std::cos(2 * M_PI * i);
//        std::array<double, 3> m = {x, y, 0};
//        return Ku1 *
//}

TEST(CubicAnisotropyField, Energy_analytic_check) {
    double Ku1 = 1e6;
    double Ku2 = 0;
    double Ku3 = 0;
    CubicAnisotropyField caniso(Ku1, Ku2, Ku3, {1, 0, 0}, {0, 1, 0});
    // Mesh mesh(1, 1, 1, 1, 1, 1);
    Mesh mesh(1, 1, 1, 1e-9, 1e-9, 1e-9);
    // double Ms = 1;
    double Ms = 1 / constants::mu0;
    af::array m = af::constant(0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(mesh, Ms, m);

    // af::array h = caniso.h(state);
    // af::print("", h);
    // double E = caniso.E(state);
    // std::cout << "E=" << E << std::endl;

    // Spin m in xy-plane
    const unsigned imax = 360;
    // for (unsigned i = 0; i < imax; ++i) {
    // for (unsigned i = 0; i < 1; ++i) {
    // for (unsigned i = 45; i < 46; ++i) {
    for (unsigned i = 0; i < imax; ++i) {
        double x = std::sin(2 * M_PI * i / imax);
        double y = std::cos(2 * M_PI * i / imax);
        state.m(0, 0, 0, 0) = x;
        state.m(0, 0, 0, 1) = y;

        std::array<double, 3> m = {x, y, 0};
        // auto E_Kc1 = [] ()
        double c1m2 = std::pow(dot_product(caniso.c1, m), 2);
        double c2m2 = std::pow(dot_product(caniso.c2, m), 2);
        double c3m2 = std::pow(dot_product(caniso.c3, m), 2);

        // double c1m = dot_product(caniso.c1, m);
        // double c2m = dot_product(caniso.c2, m);
        // double c3m = dot_product(caniso.c3, m);
        // std::cout << "c1m=" << c1m << std::endl;
        // std::cout << "c2m=" << c2m << std::endl;
        // std::cout << "c3m=" << c3m << std::endl;
        // std::cout << "c1m2=" << c1m2 << std::endl;
        // std::cout << "c2m2=" << c2m2 << std::endl;
        // std::cout << "c3m2=" << c3m2 << std::endl;

        double E_Kc1 = caniso.Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2) * (mesh.dx * mesh.dy * mesh.dz);
        // double E_Kc1 = caniso.Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2) * mesh.dx * mesh.dy * mesh.dz;
        // std::cout << E_Kc1 << ". " << caniso.E(state) / (mesh.dx * mesh.dy * mesh.dz) << std::endl;
        std::cout.precision(16);
        // std::cout << i << " " << x << " " << y << "; E:" << caniso.E(state) << " " << E_Kc1 << std::endl;
        // std::cout << i << " " << caniso.E(state) / (mesh.dx * mesh.dy * mesh.dz) << " " << E_Kc1 << std::endl;
        // std::cout << i << " " << caniso.E(state) << " " << E_Kc1 << std::endl << std::endl;
        // std::cout << i << " " << std::fabs(caniso.E(state)) * 2 << " " << E_Kc1 << std::endl;
        // EXPECT_EQ(E_Kc1, caniso.E(state));
        // EXPECT_NEAR(E_Kc1, caniso.E(state), 1e-7);
        EXPECT_NEAR(E_Kc1, caniso.E(state), 1e-37);
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
