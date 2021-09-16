// RKKY example from https://mumax.github.io/examples.html
#include "field_terms/micro/demag_field.hpp"
#include "field_terms/micro/rkky_exchange_field.hpp"
#include "integrators/llg_integrator.hpp"
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

using namespace magnumaf;

TEST(RKKY, mumax3test) {

    // Parameter initialization
    const int nx = 10, ny = 10, nz = 2;
    const double dx = 1e-9;

    const double Ms = 1e6;
    const double A = 10e-12;
    const double RKKY = -1e-3 * dx;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);

    // Initial magnetic field
    af::array m = af::constant(0.0, mesh::dims_v(mesh), f64);
    m(af::span, af::span, af::span, 0) = 1.;
    State state(mesh, Ms, m);
    af::array rkkyvals = af::constant(RKKY / 2., mesh::dims_v(mesh), f64);
    af::array exchvals = af::constant(A, mesh::dims_v(mesh), f64);
    auto rkky = RKKYExchangeField(RKKY_values(rkkyvals), Exchange_values(exchvals), mesh, false);

    auto demag = DemagField(mesh, false, false, 0);
    LLGIntegrator llg(1.0, fieldterm::to_vec(std::move(demag), std::move(rkky)));

    std::vector<double> vecE;

    for (int i = 0; i < 360; i++) {
        const double mix = std::cos(i * M_PI / 180.);
        const double miy = std::sin(i * M_PI / 180.);
        state.m(af::span, af::span, 1, 0) = mix;
        state.m(af::span, af::span, 1, 1) = miy;

        double E = llg.E(state);
        vecE.push_back(E);
        // std::cout << "i = " << i <<  ", E= " << E << std::endl;
    }

    double maxval = *std::max_element(std::begin(vecE), std::end(vecE));
    double minval = *std::min_element(std::begin(vecE), std::end(vecE));

    // std::cout << "maxval = " << maxval << std::endl;
    // std::cout << "minval = " << minval << std::endl;
    // std::cout << "Diff maxval minval  = " << maxval - minval << std::endl;
    EXPECT_NEAR(maxval - minval, 2.13934e-19, 5e-26);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
