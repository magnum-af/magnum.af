#include "mesh.hpp"
#include "micro_exch_sparse.hpp"
#include <gtest/gtest.h>

TEST(micro_exch_sparse, COO_CSR_comparison) {
    unsigned nx = 10, ny = 10, nz = 10;
    const double A = 1e-12;
    magnumafcpp::Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
    auto exch_COO = magnumafcpp::SparseExchangeField(
        af::constant(A, nx, ny, nz, f64), mesh, true, true);
    auto exch_CSR = magnumafcpp::SparseExchangeField(
        af::constant(A, nx, ny, nz, f64), mesh, true, false);
    af::array m =
        af::mean(af::matmul(exch_COO.matr - exch_CSR.matr,
                            af::flat(af::constant(1., nx, ny, nz, 3, f64))));
    // af::print("", m);
    ASSERT_EQ(m.scalar<double>(), 0);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
