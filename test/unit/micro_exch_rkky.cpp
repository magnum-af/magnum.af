#include "micro/rkky_exchange_field.hpp"
#include "util/version.hpp"
#include <gtest/gtest.h>

namespace magnumafcpp {

TEST(micro_exch_sparse, COO_CSR_comparison) {
    if (Version(af_version_string()) >= Version("3.7.0")) {
        unsigned nx = 10, ny = 10, nz = 10;
        const double A = 1e-12;
        const double RKKY = 1e-3 * 1e-9;
        Mesh mesh(nx, ny, nz, 1e-9, 1e-9, 1e-9);
        af::array rkkyvals = af::constant(RKKY / 2., mesh.dims, f64);
        af::array exchvals = af::constant(A, mesh.dims, f64);
        auto COO = RKKYExchangeField(RKKY_values(rkkyvals), Exchange_values(exchvals), mesh, af::array(), false, true);
        auto CSR = RKKYExchangeField(RKKY_values(rkkyvals), Exchange_values(exchvals), mesh, af::array(), false, false);

        af::array m = af::mean(af::matmul(COO.matr - CSR.matr, af::flat(af::constant(1., nx, ny, nz, 3, f64))));
        EXPECT_EQ(m.scalar<double>(), 0);
    }
}
} // namespace magnumafcpp

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
