#include "magnum_af.hpp"
#include "util/arg_parser.hpp"
#include <filesystem>

using namespace magnumafcpp;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;

    // Parameter initialization
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 100, ny = 25, nz = 1;

    const double A = 1.3e-11;
    const double Ms = 8e5;

    // Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);

    // Initial magnetic field
    af::array m = af::constant(0, nx, ny, nz, 3, f64);
    m(0, af::span, af::span, 1) = 1;
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(-1, af::span, af::span, 1) = 1;

    // State object
    State state(mesh, Ms, m);
    state.write_vti(outdir / "minit");

    DemagField dmag(mesh, true, true, 0);
    ExchangeField exch(A);
    LLGIntegrator Llg(1, fieldterm::to_vec(std::move(dmag), std::move(exch)));

    std::ofstream stream(outdir / "m.dat");
    stream.precision(12);

    // Relax
    StageTimer timer;
    while (state.t < 1e-9) {
        Llg.step(state);
        stream << state << std::endl;
    }
    timer.print_stage("relax ");
    state.write_vti(outdir / "relax");

    // Prepare switch
    af::array zeeswitch = af::constant(0.0, nx, ny, nz, 3, f64);
    zeeswitch(af::span, af::span, af::span, 0) = -24.6e-3 / constants::mu0;
    zeeswitch(af::span, af::span, af::span, 1) = +4.3e-3 / constants::mu0;
    Llg.llgterms.push_back(fieldterm::to_uptr<ExternalField>(zeeswitch));
    Llg.alpha = 0.02;

    // Switch
    while (state.t < 2e-9) {
        Llg.step(state);
        stream << state << std::endl;
    }
    state.write_vti(outdir / "2ns");
    stream.close();
    timer.print_stage("switch");
    timer.print_accumulated();
    return 0;
}
