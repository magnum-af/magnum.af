#include "mesh.hpp"
#include "micro_exch_sparse.hpp"
#include <gtest/gtest.h>

namespace magnumafcpp {

TEST(micro_exch_sparse, scalar_A_COO_CSR_comparison) {
    unsigned nx = 10, ny = 10, nz = 10;
    const double A = 1e-12;
    Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
    auto exch_COO = SparseExchangeField(A, mesh, true, true);
    auto exch_CSR = SparseExchangeField(A, mesh, true, false);
    af::array m =
        af::mean(af::matmul(exch_COO.matr - exch_CSR.matr,
                            af::flat(af::constant(1., nx, ny, nz, 3, f64))));
    // af::print("", m);
    ASSERT_EQ(m.scalar<double>(), 0);
}

TEST(micro_exch_sparse, array_A_COO_CSR_comparison) {
    unsigned nx = 10, ny = 10, nz = 10;
    const double A = 1e-12;
    Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
    auto exch_COO =
        SparseExchangeField(af::constant(A, nx, ny, nz, f64), mesh, true, true);
    auto exch_CSR = SparseExchangeField(af::constant(A, nx, ny, nz, f64), mesh,
                                        true, false);
    af::array m =
        af::mean(af::matmul(exch_COO.matr - exch_CSR.matr,
                            af::flat(af::constant(1., nx, ny, nz, 3, f64))));
    // af::print("", m);
    ASSERT_EQ(m.scalar<double>(), 0);
}
} // namespace magnumafcpp

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
