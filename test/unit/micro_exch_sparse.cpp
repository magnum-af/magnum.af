#include "mesh.hpp"
#include "micro/sparse_exchange_field.hpp"
#include "util/af_overloads.hpp" // TODO DEL
#include <gtest/gtest.h>

namespace magnumafcpp {

TEST(micro_exch_sparse, scalar_A_COO_CSR_comparison) {
    unsigned nx = 10, ny = 10, nz = 10;
    const double A = 1e-12;
    Mesh mesh(nx, ny, nz, 1e-9, 2e-9, 1e-9);
    auto exch_COO = SparseExchangeField(A, mesh, true, true);
    auto exch_CSR = SparseExchangeField(A, mesh, true, false);
    af::array diff =
        af::sum(af::matmul(exch_COO.get_matr() - exch_CSR.get_matr(), af::flat(af::constant(1., nx, ny, nz, 3, f64))));
    // af::print("diff", diff);
    EXPECT_EQ(diff.scalar<double>(), 0);
}

TEST(micro_exch_sparse, array_A_COO_CSR_comparison) {
    unsigned nx = 8, ny = 10, nz = 12;
    Mesh mesh(nx, ny, nz, 1e-9, 2e-9, 3e-9);
    const double A = 1e-12;
    af::array A_arr = A * af::randn(mesh::dims_s(mesh), f64);
    // af::array A_arr = af::constant(A, nx, ny, nz, f64);
    // af::print("", A_arr);
    A_arr(nx / 2, ny / 2, nz / 2) = 0.0;
    auto exch_COO = SparseExchangeField(A_arr, mesh, true, true);
    auto exch_CSR = SparseExchangeField(A_arr, mesh, true, false);
    af::array diff =
        af::sum(af::matmul(exch_COO.get_matr() - exch_CSR.get_matr(), af::flat(af::constant(1., nx, ny, nz, 3, f64))));
    EXPECT_EQ(diff.scalar<double>(), 0);
}

TEST(micro_exch_sparse, H_field_homogenuous_cube) {
    unsigned nx = 10, ny = 10, nz = 10;
    const double A = 1e-12;
    Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
    auto m = af::constant(0.0, mesh::dims_v(mesh), f64);
    m(af::span, af::span, af::span, 2) = 1.0;
    State state(mesh, 8e5, m);
    auto exch_COO = SparseExchangeField(af::constant(A, nx, ny, nz, f64), mesh, true, true);
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
    auto sum3d = [](const af::array& a) { return af::sum(af::sum(af::sum(a, 0), 1), 2).as(f64).scalar<double>(); };
    auto sum4d = [](const af::array& a) {
        return af::sum(af::sum(af::sum(af::sum(a, 0), 1), 2), 3).as(f64).scalar<double>();
    };

    EXPECT_EQ(sum3d(h(af::span, af::span, af::span, 2)), 0);
    EXPECT_EQ(sum4d(h), 0);

    state.m(af::span, af::span, af::seq(0, nz / 2 - 1), 2) = -1.0;
    auto h_m2 = exch_COO.H_in_Apm(state);
    auto HxHysum_m2 = mean_all_xy(h_m2);
    EXPECT_EQ(HxHysum_m2.scalar<double>(), 0);

    // af::print("h_m2x", h_m2(nx / 2, ny / 2, af::span, 0)); // should be zero
    // af::print("h_m2y", h_m2(nx / 2, ny / 2, af::span, 1)); // should be lero

    // af::print("h_m2", h_m2(nx / 2, ny / 2, af::span, 2)); // should zero except for boundary

    // Should be zero except for boundary at nz/2
    EXPECT_EQ(sum3d(h_m2(af::span, af::span, af::seq(0, nz / 2 - 2), 2)), 0);
    EXPECT_EQ(sum3d(h_m2(af::span, af::span, af::seq(nz / 2 + 1, nz - 1), 2)), 0);
    // Boundary
    // TODO: test exact value, not only that != zero
    EXPECT_NE(sum3d(af::abs(h_m2(af::span, af::span, af::seq(nz / 2 - 1, nz / 2), 2))), 0);
}

TEST(micro_exch_sparse, compare_a_fields) {
    // TODO make random m field and compare exch vs sparse and PBC
}

} // namespace magnumafcpp

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
