#include "field_terms/micro_exch.hpp"
#include <gtest/gtest.h>
#include <iostream>

using namespace magnumafcpp;

// Testing whether material.A and state.micro_A_field yield same result
TEST(ExchangeField, A_scalar_vs_array_value) {
    double A = 1.3e-11;
    double Ms = 8e5;
    Mesh mesh(3, 3, 3, 0.1, 0.2, 0.3);
    State state(mesh, Ms, mesh.init_sp4());
    ExchangeField exch_global(A);
    af::array globalA = exch_global.h(state);

    af::array A_field = af::constant(A, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    ExchangeField exch_local(A_field);
    af::array localA = exch_local.h(state);

    for (unsigned n0 = 0; n0 < mesh.n0; n0++) {
        for (unsigned n1 = 0; n1 < mesh.n1; n1++) {
            for (unsigned n2 = 0; n2 < mesh.n2; n2++) {
                for (unsigned n3 = 0; n3 < 3; n3++) {
                    // std::cout << afvalue(constantA(n0, n1, n2, n3)) << "\t"
                    // <<  afvalue(variableA(n0, n1, n2, n3)) << std::endl;
                    EXPECT_EQ(globalA(n0, n1, n2, n3).scalar<double>(), localA(n0, n1, n2, n3).scalar<double>());
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
