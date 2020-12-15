#include "field_terms/micro/demag_field_pbc.hpp"
#include "util/vtk_IO.hpp"
#include <chrono>
#include <gtest/gtest.h>

using namespace magnumafcpp;
auto test(Mesh mesh, double Ms = 8e5, std::size_t m_layer = 7, double opclerr = 4e-3, double err = 0.0) {
    DemagFieldPBC demag_pbc;
    const std::size_t nz = 2;
    af::dtype type = f32; // works for f64, f32, f16
    af::array m = af::constant(0, mesh::dims_v(mesh), type);
    m(af::span, af::span, m_layer, nz) = 1.;
    State state(mesh, Ms, m);
    auto H_in_Apm = demag_pbc.H_in_Apm(state);
    auto Hz_layer = H_in_Apm(mesh.nx / 2, mesh.ny / 2, m_layer, nz).as(f64).scalar<double>();
    for (std::size_t i = 0; i < mesh.nz; ++i) {
        auto Hzi = H_in_Apm(mesh.nx / 2, mesh.ny / 2, 0, nz).as(f64).scalar<double>();
        if (i != m_layer) {
            if (af::getActiveBackend() == AF_BACKEND_OPENCL) {
                EXPECT_NEAR(std::abs(Hz_layer - Hzi), Ms, opclerr);
            } else {
                // Note: cpu is more precise for f32 here, could be coincidence
                EXPECT_NEAR(std::abs(Hz_layer - Hzi), Ms, err);
            }
        }
    }
    return Hz_layer;
}

// Testing H_in_Apm of an infinite magnetic plane in xy
TEST(DemagFieldPBC, magnetic_xy_plane) {

    auto amp1 = test(Mesh{32, 18, 16, 1e-9, 1e-9, 1e-9}, 8e5, 7, 0.2); // opclerr = 4e-3 for f64
    auto amp2 = test(Mesh{32, 18, 16, 1e-9, 1e-9, 2e-9}, 8e5, 7, 4e-1, 2e-3);
    EXPECT_NEAR(amp1, amp2, 0); // TODO does not scale with 1/dz
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
