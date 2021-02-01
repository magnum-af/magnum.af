#include "field_terms/micro/demag_field_pbc.hpp"
#include <gtest/gtest.h>

using namespace magnumafcpp;
auto test(Mesh mesh, double Ms = 8e5, std::size_t m_layer = 7, std::size_t nm = 0, double gpuerr = 4e-3,
          double err = 0.0) {
    DemagFieldPBC demag_pbc;
    const std::size_t nz = 2;
    af::dtype type = f64; // works for f64, f32, f16
    af::array m = af::constant(0, mesh::dims_v(mesh), type);
    m(af::span, af::span, af::seq(m_layer, m_layer + nm), nz) = 1.;
    const State state(mesh, Ms, m);
    const auto H_in_Apm = demag_pbc.H_in_Apm(state);
    const auto Hz_layer = H_in_Apm(mesh.nx / 2, mesh.ny / 2, m_layer, nz).as(f64).scalar<double>();
    const auto Hz0 = H_in_Apm(mesh.nx / 2, mesh.ny / 2, 0, nz).as(f64).scalar<double>();

    for (std::size_t i = 0; i < mesh.nz; ++i) {
        auto Hzi = H_in_Apm(mesh.nx / 2, mesh.ny / 2, i, nz).as(f64).scalar<double>();
        // std::cout << "Hzi " << Hzi << " Hz_layer " << Hz_layer << std::endl;
        // std::cout << "Hzi " << Hzi << " Hz_layer " << Hz_layer << " " << Hz_layer / Hzi << std::endl;
        if (i < m_layer or i > m_layer + nm) {
            if (af::getActiveBackend() == AF_BACKEND_CPU) {
                EXPECT_NEAR(std::abs(Hz_layer - Hzi), Ms, err);
            } else {
                // Note: cpu is more precise for f32 here, could be coincidence
                EXPECT_NEAR(std::abs(Hz_layer - Hzi), Ms, gpuerr);
            }
        }
    }
    return std::array<double, 2>{Hz_layer, Hz0};
}

// Testing H_in_Apm of an infinite magnetic plane in xy
TEST(DemagFieldPBC, magnetic_xy_plane) {
    {
        auto [Hzl_1, Hz0_1] = test(Mesh{32, 18, 16, 1e-9, 1e-9, 1e-9}, 8e5, 7, 0, 0.2, 0.04); // opclerr = 4e-3 for f64
        auto [Hzl_2, Hz0_2] = test(Mesh{32, 18, 16, 1e-9, 1e-9, 2e-9}, 8e5, 7, 0, 4e-1, 0.04);
        EXPECT_NEAR(Hzl_1, Hzl_2, 0);
        EXPECT_NEAR(Hz0_1, Hz0_2, 0);
    }

    // 1/n scaling:
    // hzl_11 / hzl_2 == #nonmag layers 1 / #nonmag layer 2
    const auto nz = 16;
    auto [hzl_1, hz0_1] = test(Mesh{32, 18, nz, 1e-9, 1e-9, 1e-9}, 8e5, 7, 0, 0.2, 0.04); // opclerr = 4e-3 for f64
    auto [hzl_2, hz0_2] = test(Mesh{32, 18, nz, 1e-9, 1e-9, 1e-9}, 8e5, 7, 1, 0.2, 7e-2); // opclerr = 4e-3 for f64
    auto [hzl_3, hz0_3] = test(Mesh{32, 18, nz, 1e-9, 1e-9, 1e-9}, 8e5, 7, 2, 0.3, 0.2);  // opclerr = 4e-3 for f64
    auto [hzl_4, hz0_4] = test(Mesh{32, 18, nz, 1e-9, 1e-9, 1e-9}, 8e5, 7, 3, 0.2, 0.2);  // opclerr = 4e-3 for f64

    EXPECT_NEAR(hzl_1 / hzl_2, (nz - 1.) / (nz - 2.), 1e-7);
    EXPECT_NEAR(hzl_1 / hzl_3, (nz - 1.) / (nz - 3.), 2e-7);
    EXPECT_NEAR(hzl_1 / hzl_4, (nz - 1.) / (nz - 4.), 3e-7);

    // Hz_l / Hz0 = #nonmag layer/ #mag layers
    EXPECT_NEAR(std::abs(hzl_1 / hz0_1), (nz - 1.) / 1., 2e-5); // 1e-5 w.o. cuda
    EXPECT_NEAR(std::abs(hzl_2 / hz0_2), (nz - 2.) / 2., 1e-5);
    EXPECT_NEAR(std::abs(hzl_3 / hz0_3), (nz - 3.) / 3., 1e-5);
    EXPECT_NEAR(std::abs(hzl_4 / hz0_4), (nz - 4.) / 4., 1e-5);

    // std::cout << (nz - 1) / 1. << " " << std::abs(hzl_1 / hz0_1) << std::endl;
    // std::cout << (nz - 2) / 2. << " " << std::abs(hzl_2 / hz0_2) << std::endl;
    // std::cout << (nz - 3) / 3. << " " << std::abs(hzl_3 / hz0_3) << std::endl;
    // std::cout << (nz - 4) / 4. << " " << std::abs(hzl_4 / hz0_4) << std::endl;
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
