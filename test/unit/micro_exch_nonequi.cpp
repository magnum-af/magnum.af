#include "mesh.hpp"
#include "nonequi/nonequi_exchange_field.hpp"
#include "util/version.hpp"
#include <gtest/gtest.h>

namespace magnumafcpp {

TEST(micro_exch_sparse, COO_CSR_array_A) {
    if (Version(af_version_string()) >= Version("3.7.0")) {
        unsigned nx = 10, ny = 10;
        const double A = 1e-12;
        NonequiMesh mesh(nx, ny, 1e-9, 1e-9, {1e-9, 1e-9, 1e-9, 1e-9});
        const unsigned nz = mesh.nz;
        auto exch_COO = NonequiExchangeField(mesh, af::constant(A, nx, ny, nz, f64), false, true);
        auto exch_CSR = NonequiExchangeField(mesh, af::constant(A, nx, ny, nz, f64), false, false);
        af::array m =
            af::mean(af::matmul(exch_COO.matr - exch_CSR.matr, af::flat(af::constant(1., nx, ny, nz, 3, f64))));
        // af::print("", m);
        EXPECT_EQ(m.scalar<double>(), 0);
    }
}

TEST(micro_exch_sparse, COO_CSR_const_A) {
    if (Version(af_version_string()) >= Version("3.7.0")) {
        unsigned nx = 10, ny = 10;
        const double A = 1e-12;
        NonequiMesh mesh(nx, ny, 1e-9, 1e-9, {1e-9, 1e-9, 1e-9, 1e-9});
        const unsigned nz = mesh.nz;
        auto exch_COO = NonequiExchangeField(mesh, A, false, true);
        auto exch_CSR = NonequiExchangeField(mesh, A, false, false);
        af::array m =
            af::mean(af::matmul(exch_COO.matr - exch_CSR.matr, af::flat(af::constant(1., nx, ny, nz, 3, f64))));
        // af::print("", m);
        EXPECT_EQ(m.scalar<double>(), 0);
    }
}
} // namespace magnumafcpp

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
