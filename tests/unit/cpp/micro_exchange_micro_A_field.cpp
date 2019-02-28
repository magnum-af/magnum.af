#include <gtest/gtest.h>
#include <iostream>
#include "../../../src/mesh.cpp"
#include "../../../src/material.cpp"
#include "../../../src/state.cpp"
#include "../../../src/func.cpp"
#include "../../../src/vtk_IO.cpp"
#include "../../../src/misc.cpp"
#include "../../../src/integrators/adaptive_runge_kutta.cpp"
#include "../../../src/integrators/new_llg.cpp"
#include "../../../src/integrators/controller.cpp"
#include "../../../src/llg_terms/micro_exch.cpp"

// Testing whether material.A and state.micro_A_field yield same result
TEST(MicroAFieldTest, n) {
    Mesh mesh(3,3,3,0.1,0.2,0.3);
    Material material = Material();
    material.ms    = 8e5;
    material.A     = 1.3e-11;
    material.alpha = 1;
    State state(mesh,material, mesh.init_sp4());
    LlgTerms llgterm;
    llgterm.push_back(  std::shared_ptr<LLGTerm> (new ExchangeField(mesh,material)));
    LLGIntegrator Llg(llgterm);
    af::array constantA = Llg.llgterms[0]->h(state);

    state.micro_A_field=af::constant(material.A, mesh.n0, mesh.n1, mesh.n2, 3, f64);

    af::array variableA = Llg.llgterms[0]->h(state);

    for (int n0 = 0; n0 < mesh.n0; n0++){
        for (int n1 = 0; n1 < mesh.n1; n1++){
            for (int n2 = 0; n2 < mesh.n2; n2++){
                for (int n3 = 0; n3 < 3; n3++){
                    //std::cout << afvalue(constantA(n0, n1, n2, n3)) << "\t" <<  afvalue(variableA(n0, n1, n2, n3)) << std::endl;
                    ASSERT_EQ(afvalue(constantA(n0, n1, n2, n3)), afvalue(variableA(n0, n1, n2, n3)));
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
