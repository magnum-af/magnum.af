// ref(section 4., Fig 1.): D. Corts-Ortuo et al. , “Proposal for a
// micromagnetic standard problem for materials with Dzyaloshinskii-Moriya
// interaction,” New Journal of Physics, vol. 20, p. 113015, Nov. 2018.
#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    // Checking input variables and setting GPU Device
    for (int i = 0; i < argc; i++) {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();

    // Parameter initialization
    const double x = 1.e-7, y = 1e-9, z = 1e-9;
    const int nx = 100, ny = 1, nz = 1;
    const double A = 13e-12;
    const double D = 3e-3;
    const double Ms = 0.86e6;
    const double Ku = 0.4e6;

    // Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);

    // Initial magnetic field
    af::array m = af::constant(0, mesh.dims, f64);
    m(af::span, af::span, af::span, 2) = 1;
    State state(mesh, Ms, m);
    state.write_vti(filepath + "minit");

    // auto demag = LlgTerm (new DemagField(mesh, true, true, 0));
    auto exch = LlgTerm(new ExchangeField(A));
    auto aniso = LlgTerm(new UniaxialAnisotropyField(Ku, {0, 0, 1}));
    auto dmi = LlgTerm(new DmiField(D, {0, 0, 1}));
    LLGIntegrator Llg(1, {exch, aniso, dmi});

    std::ofstream stream;
    stream.precision(12);
    stream.open(filepath + "m.dat");

    // Relax
    StageTimer timer;
    while (state.t < 2e-10) {
        Llg.step(state);
        stream << state << std::endl;
        if (state.steps % 100 == 0)
            state.write_vti(filepath + "m_step" + std::to_string(state.steps));
    }
    // Llg.relax(state, 1e-10, 100, 1);
    stream.close();
    timer.print_stage("relax ");
    state.write_vti(filepath + "relax");

    // write m to file
    // plot with; gnuplot -e 'set terminal pdf; set output "m.pdf"; set xlabel
    // "x [nm]"; set ylabel "m_i"; p "mrelaxed.dat" u 1:2 w lp t "mx", "" u 1:4
    // w lp t "mz"'
    std::ofstream mstream(filepath + "mrelaxed.dat");
    mstream << "# idx, mx, my, mz" << std::endl;
    for (int i = 0; i < nx; i++) {
        mstream << i << "\t" << state.m(i, 0, 0, 0).scalar<double>() << "\t" << state.m(i, 0, 0, 1).scalar<double>()
                << "\t" << state.m(i, 0, 0, 2).scalar<double>() << std::endl;
    }
    mstream.close();

    timer.print_accumulated();
    return 0;
}
