#include "../../src/mesh.hpp"
#include "../../src/nonequispaced_mesh.hpp"//TODO needed for vtk_IO, include headers instead?
#include "../../src/material.hpp"
#include "../../src/state.hpp"
#include "../../src/func.hpp"
#include "../../src/vtk_IO.hpp"
#include "../../src/misc.hpp"
#include "../../src/llg_terms/micro_exch.hpp"
#include <gtest/gtest.h>
#include <iostream>

using namespace magnumafcpp;


// Testing whether material.A and state.micro_A_field yield same result
TEST(StateMicroAField, MicroASingleValueVsArrayHeffTest) {
    double A = 1.3e-11;
    double Ms    = 8e5;
    Mesh mesh(3, 3, 3, 0.1, 0.2, 0.3);
    State state(mesh, Ms, mesh.init_sp4());
    ExchangeField exch_global(A);
    af::array globalA = exch_global.h(state);

    af::array A_field = af::constant(A, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    ExchangeField exch_local(A_field);
    af::array localA = exch_local.h(state);

    for (uint32_t n0 = 0; n0 < mesh.n0; n0++){
        for (uint32_t n1 = 0; n1 < mesh.n1; n1++){
            for (uint32_t n2 = 0; n2 < mesh.n2; n2++){
                for (uint32_t n3 = 0; n3 < 3; n3++){
                    //std::cout << afvalue(constantA(n0, n1, n2, n3)) << "\t" <<  afvalue(variableA(n0, n1, n2, n3)) << std::endl;
                    ASSERT_EQ(globalA(n0, n1, n2, n3).scalar<double>(), localA(n0, n1, n2, n3).scalar<double>());
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
