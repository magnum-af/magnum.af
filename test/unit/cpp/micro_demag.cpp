#include "../../../src/mesh.cpp"
#include "../../../src/nonequispaced_mesh.cpp"//TODO needed for vtk_IO, include headers instead?
#include "../../../src/llg_terms/micro_demag.cpp"
#include "../../../src/integrators/new_llg.cpp"
#include "../../../src/integrators/adaptive_runge_kutta.cpp"
#include "../../../src/integrators/controller.cpp"
#include "../../../src/func.cpp"
#include "../../../src/state.cpp"
#include "../../../src/misc.cpp"
#include "../../../src/vtk_IO.cpp"
#include <gtest/gtest.h>
using namespace magnumafcpp;

TEST(MicroDemag, EnergyOfHomogeneousCube) {
    const double x=1.e-9, y=1.e-9, z=1.e-9;
    const int nx = 10, ny=10 , nz=10;

    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    double Ms = 1e5;

    af::array m = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(af::span, af::span, af::span, 0) = 1;
    State state(mesh, Ms, m);

    LLGIntegrator llg(1., {LlgTerm (new DemagField(mesh))});

    std::cout.precision(24);
    double llgE = llg.E(state);
    double analytic = 1./6. * x * y * z * pow(state.Ms, 2) * constants::mu0;

    EXPECT_NEAR(llgE, analytic, (llgE+analytic)/2. * 5e-8);// opencl precision
    //EXPECT_NEAR(llgE, analytic, (llgE+analytic)/2. * 1e-12);// cpu precision
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
