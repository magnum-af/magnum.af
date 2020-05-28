#include "micro_exch_nonequi.hpp"
#include "mesh.hpp"
#include <gtest/gtest.h>

namespace magnumafcpp {

TEST(micro_exch_sparse, COO_CSR_array_A) {
    unsigned nx = 10, ny = 10;
    const double A = 1e-12;
    NonequispacedMesh mesh(nx, ny, 1e-9, 1e-9, {1e-9, 1e-9, 1e-9, 1e-9});
    const unsigned nz = mesh.nz;
    auto exch_COO = NonequiExchangeField(af::constant(A, nx, ny, nz, f64), mesh,
                                         false, true);
    auto exch_CSR = NonequiExchangeField(af::constant(A, nx, ny, nz, f64), mesh,
                                         false, false);
    af::array m =
        af::mean(af::matmul(exch_COO.matr - exch_CSR.matr,
                            af::flat(af::constant(1., nx, ny, nz, 3, f64))));
    // af::print("", m);
    ASSERT_EQ(m.scalar<double>(), 0);
}

TEST(micro_exch_sparse, COO_CSR_const_A) {
    unsigned nx = 10, ny = 10;
    const double A = 1e-12;
    NonequispacedMesh mesh(nx, ny, 1e-9, 1e-9, {1e-9, 1e-9, 1e-9, 1e-9});
    const unsigned nz = mesh.nz;
    auto exch_COO = NonequiExchangeField(A, mesh, false, true);
    auto exch_CSR = NonequiExchangeField(A, mesh, false, false);
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
