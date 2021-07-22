#include "mesh.hpp"
#include "micro/exchange_field_pbc.hpp"
#include <gtest/gtest.h>

namespace magnumafcpp {

TEST(micro_exch_pbc, scalar_A_COO_CSR_comparison) {
    unsigned nx = 10, ny = 10, nz = 10;
    const double A = 1e-12;
    Mesh mesh(nx, ny, nz, 1e-9, 2e-9, 1e-9);
    auto exch_COO = ExchangeFieldPBC(A, mesh, true, true);
    auto exch_CSR = ExchangeFieldPBC(A, mesh, true, false);
    af::array diff =
        af::sum(af::matmul(exch_COO.get_matr() - exch_CSR.get_matr(), af::flat(af::constant(1., nx, ny, nz, 3, f64))));
    // af::print("diff", diff);
    EXPECT_EQ(diff.scalar<double>(), 0);
}

TEST(micro_exch_pbc, array_A_COO_CSR_comparison) {
        unsigned nx = 10, ny = 10, nz = 10;
        const double A = 1e-12;
        Mesh mesh(nx, ny, nz, 1e-9, 2e-9, 3e-9);
        af::array A_arr = af::constant(A, nx, ny, nz, f64);
        A_arr(nx / 2, ny / 2, nz / 2) = 0.0;
        auto exch_COO = ExchangeFieldPBC(A_arr, mesh, true, true);
        auto exch_CSR = ExchangeFieldPBC(A_arr, mesh, true, false);
        af::array diff = af::sum(
            af::matmul(exch_COO.get_matr() - exch_CSR.get_matr(), af::flat(af::constant(1., nx, ny, nz, 3, f64))));
        EXPECT_EQ(diff.scalar<double>(), 0);
}

TEST(micro_exch_pbc, H_field_homogenuous_cube) {
    unsigned nx = 10, ny = 10, nz = 10;
    const double A = 1e-12;
    Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
    auto m = af::constant(0.0, mesh::dims_v(mesh), f64);
    m(af::span, af::span, af::span, 2) = 1.0;
    State state(mesh, 8e5, m);
    auto exch_COO = ExchangeFieldPBC(af::constant(A, nx, ny, nz, f64), mesh, true, true);

    auto sum3d = [](const af::array& a) { return af::sum(af::sum(af::sum(a, 0), 1), 2).as(f64).scalar<double>(); };
    {
        auto h = exch_COO.H_in_Apm(state);
        // check if all Hx and Hy are zero
        auto mean_all_xy = [](const auto& a) {
            return af::sum(af::sum(af::sum(af::sum(a(af::span, af::span, af::span, af::seq(0, 1)), 0), 1), 2), 3);
        };
        auto HxHysum = mean_all_xy(h);
        // af::print("h", h(nx / 2, ny / 2, af::span, 0)); // should be 0
        // af::print("h", h(nx / 2, ny / 2, af::span, 1)); // should be 0
        EXPECT_EQ(HxHysum.scalar<double>(), 0);

        // With -6 * m_i, Hz should be zero as well:
        // af::print("h", h(nx / 2, ny / 2, af::span, 2)); // should be lapM
        EXPECT_EQ(sum3d(h(af::span, af::span, af::span, 2)), 0);

        // setting m half cube to -z
        state.m(af::span, af::span, af::seq(0, nz / 2 - 1), 2) = -1.0;
        auto h_m2 = exch_COO.H_in_Apm(state);
        auto HxHysum_m2 = mean_all_xy(h_m2);
        // af::print("h_m2", h_m2(nx / 2, ny / 2, af::span, 0)); // should be lap m
        // af::print("h_m2", h_m2(nx / 2, ny / 2, af::span, 1)); // should be lap m
        EXPECT_EQ(HxHysum_m2.scalar<double>(), 0);

        // H should be zero, except for the boundaries at nz-1/0 and nz/2-1  nz/2
        // af::print("h_m2", h_m2(nx / 2, ny / 2, af::span, 2)); // should be lap m
        EXPECT_EQ(sum3d(h_m2(af::span, af::span, af::seq(1, nz / 2 - 2), 2)), 0);
        EXPECT_EQ(sum3d(h_m2(af::span, af::span, af::seq(nz / 2 + 1, nz - 2), 2)), 0);
        // Looks good, affects over boundary:
        // TODO check amplitude vs analytical (for all exch impls)
        // TODO check exact numeric values
        // TODO check amplitude at 0, z/2 and z boundaries
        EXPECT_NE(sum3d(af::abs(h_m2(af::span, af::span, 0, 2))), 0);
        EXPECT_NE(sum3d(af::abs(h_m2(af::span, af::span, af::seq(nz / 2 - 1, nz / 2), 2))), 0);
        EXPECT_NE(sum3d(af::abs(h_m2(af::span, af::span, nz - 1, 2))), 0);


    }
    {
        // Looks good, affects over boundary:
        state.m = 0;
        state.m(af::span, af::span, 0, 2) = 1.0;
        auto h_m3 = exch_COO.H_in_Apm(state);
        af::print("h_m3", h_m3(nx / 2, ny / 2, af::span, 2));
        EXPECT_EQ(sum3d(h_m3(af::span, af::span, af::seq(2, nz - 2), 2)), 0);
        // TODO test amplitudes
        EXPECT_NE(sum3d(af::abs(h_m3(af::span, af::span, af::seq(0, 1), 2))), 0);
        EXPECT_NE(sum3d(af::abs(h_m3(af::span, af::span, nz - 1, 2))), 0);
    }
}

// TEST(micro_exch_pbc, timing) {
//     unsigned nx = 100, ny = 100, nz = 100;
//     const double A = 1e-12;
//     Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
//     auto exch_COO = ExchangeFieldPBC(af::constant(A, nx, ny, nz, f64), mesh, true, true);
//     auto exch_CSR = ExchangeFieldPBC(af::constant(A, nx, ny, nz, f64), mesh, true, false);
//
//     auto m = af::constant(0.0, mesh::dims_v(mesh), f64);
//     m(af::span, af::span, af::span, 2) = 1.0;
//     State state(mesh, 8e5, m);
//     auto h = exch_COO.H_in_Apm(state);
// }

} // namespace magnumafcpp

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
