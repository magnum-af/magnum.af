#include "micro/cubic_anisotropy_field.hpp"
#include "util/util.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
using namespace magnumaf;
TEST(CubicAnisotropyField, Constructor) {
    CubicAnisotropyField caniso(1, 1, 1, {1, 1, 0}, {-1, 1, 0});
    const auto c1 = std::get<std::array<double, 3>>(caniso.c1.variant);
    const auto c2 = std::get<std::array<double, 3>>(caniso.c2.variant);
    const auto c3 = std::get<std::array<double, 3>>(caniso.c3.variant);

    EXPECT_THAT(c1, testing::ElementsAre(1 / std::sqrt(2), 1 / std::sqrt(2), 0));
    EXPECT_THAT(c2, testing::ElementsAre(-1 / std::sqrt(2), 1 / std::sqrt(2), 0));
    EXPECT_THAT(c3, testing::ElementsAre(0, 0, 1));
}

// single spin rotating 360 degree in xy plane with offset of z in z-dir (befor normalization)
// comparing calculated energy vs analytical for Kc1, Kc2 and Kc3 separately
void energy_test_xy_rotation(const double Kc1, const double Kc2, const double Kc3, const double z) {
    const Mesh mesh(2, 2, 2, 1e-9, 1e-9, 1e-9);
    // const Mesh mesh(1, 1, 1, 1e-9, 1e-9, 1e-9); // Alternative
    const unsigned num_of_cells = mesh.nx * mesh.ny * mesh.nz;
    // const double Ms = 1 / constants::mu0; // Alternative
    af::array Ms = af::constant(1 / constants::mu0, mesh::dims_s(mesh), f64);
    af::array m = af::constant(0, mesh::dims_v(mesh), f64);
    m(af::span, af::span, af::span, 0) = 1;
    State state(mesh, Ms, m);

    // Spin m in xy-plane
    const unsigned imax = 64;
    for (unsigned i = 0; i < imax; ++i) {
        double x = std::sin(2 * M_PI * i / imax);
        double y = std::cos(2 * M_PI * i / imax);
        state.m(af::span, af::span, af::span, 0) = x;
        state.m(af::span, af::span, af::span, 1) = y;
        state.m(af::span, af::span, af::span, 2) = z;
        util::normalize_inplace(state.m);

        const std::array<double, 3> m = util::normalize_vector({x, y, z});
        const std::array<double, 3> c1 = {1, 0, 0};
        const std::array<double, 3> c2 = {0, 1, 0};
        const std::array<double, 3> c3 = {0, 0, 1};

        const double c1m2 = std::pow(util::dot_product(c1, m), 2);
        const double c2m2 = std::pow(util::dot_product(c2, m), 2);
        const double c3m2 = std::pow(util::dot_product(c3, m), 2);
        std::cout.precision(16);

        // Kc1
        {
            CubicAnisotropyField Kc1_caniso(Kc1, 0, 0, c1, c2);
            const double Kc1_E_density_analytic = num_of_cells * Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2);
            const double Kc1_E_analytic = Kc1_E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
            const double Kc1_E_calculated = Kc1_caniso.Energy_in_J(state);
            const double Kc1_E_density_calculated = Kc1_E_calculated / (mesh.dx * mesh.dy * mesh.dz);
            EXPECT_NEAR(Kc1_E_calculated, Kc1_E_analytic, 2e-36);
            EXPECT_NEAR(Kc1_E_density_calculated, Kc1_E_density_analytic, 2e-9);
            // std::cout << i << " " << Kc1_E_density_analytic << " " << Kc1_E_density_calculated;
            // std::cout << " " << Kc1_E_analytic << " " << Kc1_E_calculated << std::endl;
        }

        // Kc1 array input
        {
            auto Kc1_ = af::constant(Kc1, mesh::dims_s(mesh), f64);
            auto Kc2_ = af::constant(0, mesh::dims_s(mesh), f64);
            auto Kc3_ = af::constant(0, mesh::dims_s(mesh), f64);
            auto c1_ = af::constant(0, mesh::dims_v(mesh), f64);
            c1_(af::span, af::span, af::span, 0) = 1;
            auto c2_ = af::constant(0, mesh::dims_v(mesh), f64);
            c2_(af::span, af::span, af::span, 1) = 1;
            CubicAnisotropyField Kc1_caniso(Kc1_, Kc2_, Kc3_, c1_, c2_);
            const double Kc1_E_density_analytic = num_of_cells * Kc1 * (c1m2 * c2m2 + c1m2 * c3m2 + c2m2 * c3m2);
            const double Kc1_E_analytic = Kc1_E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
            const double Kc1_E_calculated = Kc1_caniso.Energy_in_J(state);
            const double Kc1_E_density_calculated = Kc1_E_calculated / (mesh.dx * mesh.dy * mesh.dz);
            EXPECT_NEAR(Kc1_E_calculated, Kc1_E_analytic, 2e-36);
            EXPECT_NEAR(Kc1_E_density_calculated, Kc1_E_density_analytic, 3e-9);
            // std::cout << i << " " << Kc1_E_density_analytic << " " << Kc1_E_density_calculated;
            // std::cout << " " << Kc1_E_analytic << " " << Kc1_E_calculated << std::endl;
        }

        // Kc2, this is always zero if z == 0
        {
            CubicAnisotropyField Kc2_caniso(0, Kc2, 0, c1, c2);
            const double Kc2_E_density_analytic = num_of_cells * Kc2 * (c1m2 * c2m2 * c3m2);
            const double Kc2_E_analytic = Kc2_E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
            const double Kc2_E_calculated = Kc2_caniso.Energy_in_J(state);
            const double Kc2_E_density_calculated = Kc2_E_calculated / (mesh.dx * mesh.dy * mesh.dz);
            EXPECT_NEAR(Kc2_E_calculated, Kc2_E_analytic, 1e-36);
            EXPECT_NEAR(Kc2_E_density_calculated, Kc2_E_density_analytic, 1e-10);
            // std::cout << i << " " << Kc2_E_density_analytic << " " << Kc2_E_density_calculated;
            // std::cout << " " << Kc2_E_analytic << " " << Kc2_E_calculated << std::endl;
        }

        // Kc2 array input, this is always zero if z == 0
        {
            auto Kc1_ = af::constant(0, mesh::dims_s(mesh), f64);
            auto Kc2_ = af::constant(Kc2, mesh::dims_s(mesh), f64);
            auto Kc3_ = af::constant(0, mesh::dims_s(mesh), f64);
            auto c1_ = af::constant(0, mesh::dims_v(mesh), f64);
            c1_(af::span, af::span, af::span, 0) = 1;
            auto c2_ = af::constant(0, mesh::dims_v(mesh), f64);
            c2_(af::span, af::span, af::span, 1) = 1;
            CubicAnisotropyField Kc2_caniso(Kc1_, Kc2_, Kc3_, c1_, c2_);
            const double Kc2_E_density_analytic = num_of_cells * Kc2 * (c1m2 * c2m2 * c3m2);
            const double Kc2_E_analytic = Kc2_E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
            const double Kc2_E_calculated = Kc2_caniso.Energy_in_J(state);
            const double Kc2_E_density_calculated = Kc2_E_calculated / (mesh.dx * mesh.dy * mesh.dz);
            EXPECT_NEAR(Kc2_E_calculated, Kc2_E_analytic, 1e-36);
            EXPECT_NEAR(Kc2_E_density_calculated, Kc2_E_density_analytic, 1e-10);
            // std::cout << i << " " << Kc2_E_density_analytic << " " << Kc2_E_density_calculated;
            // std::cout << " " << Kc2_E_analytic << " " << Kc2_E_calculated << std::endl;
        }

        // Kc3
        {
            CubicAnisotropyField Kc3_caniso(0, 0, Kc3, c1, c2);
            const double c1m4 = std::pow(util::dot_product(c1, m), 4);
            const double c2m4 = std::pow(util::dot_product(c2, m), 4);
            const double c3m4 = std::pow(util::dot_product(c3, m), 4);
            const double Kc3_E_density_analytic = num_of_cells * Kc3 * (c1m4 * c2m4 + c1m4 * c3m4 + c2m4 * c3m4);
            const double Kc3_E_analytic = Kc3_E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
            const double Kc3_E_calculated = Kc3_caniso.Energy_in_J(state);
            const double Kc3_E_density_calculated = Kc3_E_calculated / (mesh.dx * mesh.dy * mesh.dz);
            EXPECT_NEAR(Kc3_E_calculated, Kc3_E_analytic, 5e-37);
            EXPECT_NEAR(Kc3_E_density_calculated, Kc3_E_density_analytic, 5e-10);
            // std::cout << i << " " << Kc3_E_density_analytic << " " << Kc3_E_density_calculated;
            // std::cout << " " << Kc3_E_analytic << " " << Kc3_E_calculated << std::endl;
        }

        // Kc3 array input
        {
            auto Kc1_ = af::constant(0, mesh::dims_s(mesh), f64);
            auto Kc2_ = af::constant(0, mesh::dims_s(mesh), f64);
            auto Kc3_ = af::constant(Kc3, mesh::dims_s(mesh), f64);
            auto c1_ = af::constant(0, mesh::dims_v(mesh), f64);
            c1_(af::span, af::span, af::span, 0) = 1;
            auto c2_ = af::constant(0, mesh::dims_v(mesh), f64);
            c2_(af::span, af::span, af::span, 1) = 1;
            CubicAnisotropyField Kc3_caniso(Kc1_, Kc2_, Kc3_, c1_, c2_);
            const double c1m4 = std::pow(util::dot_product(c1, m), 4);
            const double c2m4 = std::pow(util::dot_product(c2, m), 4);
            const double c3m4 = std::pow(util::dot_product(c3, m), 4);
            const double Kc3_E_density_analytic = num_of_cells * Kc3 * (c1m4 * c2m4 + c1m4 * c3m4 + c2m4 * c3m4);
            const double Kc3_E_analytic = Kc3_E_density_analytic * (mesh.dx * mesh.dy * mesh.dz);
            const double Kc3_E_calculated = Kc3_caniso.Energy_in_J(state);
            const double Kc3_E_density_calculated = Kc3_E_calculated / (mesh.dx * mesh.dy * mesh.dz);
            EXPECT_NEAR(Kc3_E_calculated, Kc3_E_analytic, 5e-37);
            EXPECT_NEAR(Kc3_E_density_calculated, Kc3_E_density_analytic, 5e-10);
            // std::cout << i << " " << Kc3_E_density_analytic << " " << Kc3_E_density_calculated;
            // std::cout << " " << Kc3_E_analytic << " " << Kc3_E_calculated << std::endl;
        }
    }
}

TEST(CubicAnisotropyField, energy_test_xy_rotation) {
    energy_test_xy_rotation(1e6, 1e6, 1e6, 0);
    energy_test_xy_rotation(-1e6, -1e6, -1e6, 0);
    energy_test_xy_rotation(1e6, 1e6, 1e6, 0.5);
    energy_test_xy_rotation(1e6, 1e6, 1e6, 1.0);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
