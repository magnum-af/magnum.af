#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    bool conv = false;

    for (int i = 0; i < argc; i++) {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();

    // z_y plane
    // note: m_x has negative sign compared to x_y plane
    {
        // Checking input variables and setting GPU Device
        af::timer total_time = af::timer::start();
        // Parameter initialization
        const double z = 5.e-7, y = 1.25e-7, x = 3.e-9;
        const int nx = 1, ny = 25, nz = 100;
        const double A = 1.3e-11;

        // Generating Objects
        std::vector<double> dz_spacing;
        for (int i = 0; i < nz; i++) {
            dz_spacing.push_back(z / nz);
        }
        Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
        NonequispacedMesh ne_mesh(nx, ny, x / nx, y / ny, dz_spacing);

        // Initial magnetic field
        af::array m = af::constant(0.0, nx, ny, nz, 3, f64);
        m(af::span, af::span, af::seq(1, af::end - 1), 2) = 1;
        m(af::span, af::span, 0, 1) = 1;
        m(af::span, af::span, -1, 1) = 1;

        State state(mesh, 8e5, m);
        state.nonequimesh = ne_mesh; // TODO avoid
        vti_writer_micro(state.m, mesh, (filepath + "z_minit").c_str());

        LlgTerms llgterm;
        llgterm.push_back(LlgTerm(new DemagField(mesh, true, true, 0)));
        if (conv)
            llgterm.push_back(LlgTerm(new ExchangeField(A)));
        else
            llgterm.push_back(LlgTerm(new NonequiExchangeField(ne_mesh, af::constant(A, mesh.dims, f64))));
        LLGIntegrator Llg(1, llgterm);

        std::ofstream stream;
        stream.precision(12);
        stream.open((filepath + "z_m.dat").c_str());

        // Relax
        af::timer t = af::timer::start();
        while (state.t < 1e-9) {
            Llg.step(state);
            state.calc_mean_m(stream);
        }
        std::cout << "timerelax [af-s]: " << af::timer::stop(t) << std::endl;
        vti_writer_micro(state.m, mesh, (filepath + "z_relax").c_str());

        // Prepare switch
        af::array zeeswitch = af::constant(0.0, 1, 1, 1, 3, f64);
        zeeswitch(0, 0, 0, 2) = -24.6e-3 / constants::mu0;
        zeeswitch(0, 0, 0, 1) = +4.3e-3 / constants::mu0;
        zeeswitch = tile(zeeswitch, mesh.n0, mesh.n1, mesh.n2);
        Llg.llgterms.push_back(LlgTerm(new ExternalField(zeeswitch)));
        Llg.alpha = 0.02;

        // Switch
        t = af::timer::start();
        while (state.t < 2e-9) {
            Llg.step(state);
            state.calc_mean_m(stream);
        }
        std::cout << "time integrate 1ns [af-s]: " << af::timer::stop(t) << std::endl;
        vti_writer_micro(state.m, mesh, (filepath + "2ns").c_str());
        stream.close();
        std::cout << "total [af-s]: " << af::timer::stop(total_time) << std::endl;
    }
    // x_y plane
    if (!exists(filepath + "x_m.dat")) {
        // Checking input variables and setting GPU Device
        af::timer total_time = af::timer::start();

        // Parameter initialization
        const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
        const int nx = 100, ny = 25, nz = 1;
        const double A = 1.3e-11;

        // Generating Objects
        Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);

        // Initial magnetic field
        State state(mesh, 8e5, mesh.init_sp4());
        vti_writer_micro(state.m, mesh, (filepath + "x_minit").c_str());

        LlgTerms llgterm;
        llgterm.push_back(LlgTerm(new DemagField(mesh, true, true, 0)));
        if (conv)
            llgterm.push_back(LlgTerm(new ExchangeField(A)));
        else
            llgterm.push_back(LlgTerm(new SparseExchangeField(A, mesh)));
        LLGIntegrator Llg(1, llgterm);

        std::ofstream stream;
        stream.precision(12);
        stream.open((filepath + "x_m.dat").c_str());

        // Relax
        af::timer t = af::timer::start();
        while (state.t < 1e-9) {
            Llg.step(state);
            state.calc_mean_m(stream);
        }
        std::cout << "timerelax [af-s]: " << af::timer::stop(t) << std::endl;
        vti_writer_micro(state.m, mesh, (filepath + "x_relax").c_str());

        // Prepare switch
        af::array zeeswitch = af::constant(0.0, 1, 1, 1, 3, f64);
        zeeswitch(0, 0, 0, 0) = -24.6e-3 / constants::mu0;
        zeeswitch(0, 0, 0, 1) = +4.3e-3 / constants::mu0;
        zeeswitch = tile(zeeswitch, mesh.n0, mesh.n1, mesh.n2);
        Llg.llgterms.push_back(LlgTerm(new ExternalField(zeeswitch)));
        Llg.alpha = 0.02;

        // Switch
        t = af::timer::start();
        while (state.t < 2e-9) {
            Llg.step(state);
            state.calc_mean_m(stream);
        }
        std::cout << "time integrate 1ns [af-s]: " << af::timer::stop(t) << std::endl;
        vti_writer_micro(state.m, mesh, (filepath + "2ns").c_str());
        stream.close();
        std::cout << "total [af-s]: " << af::timer::stop(total_time) << std::endl;
    }
    return 0;
}
