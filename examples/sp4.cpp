#include "magnum_af.hpp"
#include <filesystem>

using namespace magnumafcpp;

int main(int argc, char** argv) {
    // Checking input variables and setting GPU Device
    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "]= " << argv[i] << std::endl;
    }
    std::string filepath(argc > 1 ? argv[1] + std::string("/") : "output_magnum.af/");
    std::filesystem::create_directories(filepath);
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();

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
    state.write_vti(filepath + "minit");

    auto demag = uptr_FieldTerm(std::make_unique<DemagField>(mesh, true, true, 0));
    auto exch = uptr_FieldTerm(std::make_unique<ExchangeField>(A));
    LLGIntegrator Llg(1, {std::move(demag), std::move(exch)});

    std::ofstream stream;
    stream.precision(12);
    stream.open(filepath + "m.dat");

    // Relax
    StageTimer timer;
    while (state.t < 1e-9) {
        Llg.step(state);
        stream << state << std::endl;
    }
    timer.print_stage("relax ");
    state.write_vti(filepath + "relax");

    // Prepare switch
    af::array zeeswitch = af::constant(0.0, nx, ny, nz, 3, f64);
    zeeswitch(af::span, af::span, af::span, 0) = -24.6e-3 / constants::mu0;
    zeeswitch(af::span, af::span, af::span, 1) = +4.3e-3 / constants::mu0;
    Llg.llgterms.push_back(uptr_FieldTerm(new ExternalField(zeeswitch)));
    Llg.alpha = 0.02;

    // Switch
    while (state.t < 2e-9) {
        Llg.step(state);
        stream << state << std::endl;
    }
    state.write_vti(filepath + "2ns");
    stream.close();
    timer.print_stage("switch");
    timer.print_accumulated();
    return 0;
}
