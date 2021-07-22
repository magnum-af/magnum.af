#include "field_terms/micro/exchange_field.hpp"
#include "util/geometry.hpp"
#include <gtest/gtest.h>
#include <iostream>

using namespace magnumafcpp;

// Testing whether material.A and state.micro_A_field yield same result
TEST(ExchangeField, A_scalar_vs_array_value) {
    double A = 1.3e-11;
    double Ms = 8e5;
    Mesh mesh(3, 3, 3, 0.1, 0.2, 0.3);
    State state(mesh, Ms, util::init_sp4(mesh));
    ExchangeField exch_global(A);
    af::array globalA = exch_global.H_in_Apm(state);

    af::array A_field = af::constant(A, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    ExchangeField exch_local(A_field);
    af::array localA = exch_local.H_in_Apm(state);

    for (unsigned nx = 0; nx < mesh.nx; nx++) {
        for (unsigned ny = 0; ny < mesh.ny; ny++) {
            for (unsigned nz = 0; nz < mesh.nz; nz++) {
                for (unsigned n3 = 0; n3 < 3; n3++) {
                    // std::cout << util::afvalue_as_f64(constantA(nx, ny, nz, n3)) << "\t"
                    // <<  util::afvalue_as_f64(variableA(nx, ny, nz, n3)) << std::endl;
                    EXPECT_EQ(globalA(nx, ny, nz, n3).scalar<double>(), localA(nx, ny, nz, n3).scalar<double>());
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
