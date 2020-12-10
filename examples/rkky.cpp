// RKKY example from https://mumax.github.io/examples.html
#include "field_terms/micro/demag_field.hpp"
#include "field_terms/micro/exchange_field.hpp"
#include "field_terms/micro/rkky_exchange_field.hpp"
#include "integrators/llg_integrator.hpp"
#include "util/arg_parser.hpp"
#include <cmath>

using namespace magnumafcpp;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;

    // Parameter initialization
    const int nx = 10, ny = 10, nz = 2;
    const double dx = 1e-9;

    const double Ms = 1e6;
    const double A = 10e-12;
    const double RKKY = -1e-3 * dx;

    // Generating Objects
    Mesh mesh{nx, ny, nz, dx, dx, dx};

    // Initial magnetic field
    af::array m = af::constant(0.0, mesh::dims_v(mesh), f64);
    m(af::span, af::span, af::span, 0) = 1.;
    State state(mesh, Ms, m);
    state.write_vti(outdir / "minit");
    af::array rkkyvals = af::constant(RKKY / 2., mesh::dims_v(mesh), f64);
    af::array exchvals = af::constant(A, mesh::dims_v(mesh), f64);
    auto rkky = uptr_FieldTerm(new RKKYExchangeField(RKKY_values(rkkyvals), Exchange_values(exchvals), mesh));

    auto demag = uptr_FieldTerm(new DemagField(mesh, true, true, 0));
    LLGIntegrator Llg(1, {std::move(demag), std::move(rkky)});

    std::ofstream stream(outdir / "m.dat");
    stream.precision(12);

    std::ofstream streamE(outdir / "E.dat");
    streamE.precision(12);

    std::vector<double> vecE;

    for (int i = 0; i < 360; i++) {
        const double mix = std::cos(i * M_PI / 180.);
        const double miy = std::sin(i * M_PI / 180.);
        state.m(af::span, af::span, 1, 0) = mix;
        state.m(af::span, af::span, 1, 1) = miy;

        double E = Llg.E(state);
        vecE.push_back(E);
        std::cout << "i = " << i << ", E= " << E << std::endl;
        streamE << i << "\t" << E << std::endl;
        stream << state << std::endl;
        state.write_vti(outdir / ("m_field_" + std::to_string(i)));
    }
    stream.close();
    streamE.close();

    double maxval = *std::max_element(std::begin(vecE), std::end(vecE));
    double minval = *std::min_element(std::begin(vecE), std::end(vecE));

    std::cout << "maxval = " << maxval << std::endl;
    std::cout << "minval = " << minval << std::endl;
    std::cout << "Diff maxval minval  = " << maxval - minval << std::endl;

    return 0;
}
