#include "arrayfire.h"
#include "magnum_af.hpp"
#include <fstream>
#include <iostream>

using namespace magnumafcpp;

int main(int argc, char** argv) {
    std::cout << "argc" << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc > 1 ? argv[1] : "../Data/Testing");
    std::cout << "Writing into path " << filepath << std::endl;
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);

    af::info();
    std::cout.precision(24);

    std::ofstream stream(filepath + "m.dat");
    stream.precision(24);
    stream << "# nz	<mx>    <my>    <mz>    hzee    demagx  demagy demagz" << std::endl;
    // Parameter initialization
    for (int nz = 2; nz < 22; nz++) {
        std::cout << "starting nz = " << nz << std::endl;
        const int nx = 250, ny = 250; //, nz = 1;
        const double x = 1000e-9,
                     y = 1000e-9; //, z = 65e-9; //[m] // Physical dimensions

        // Generating Objects
        Mesh mesh(nx, ny, nz, x / nx, y / ny, 5e-9);
        double Ms = 1.393e6; //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
        double A = 1.5e-11;  //[J/m]
        double zee = 40e-3 / constants::mu0;

        af::array m0 = mesh.init_vortex();
        m0(af::span, af::span, 0, af::span) = 0;
        State state(mesh, Ms, m0);
        state.write_vti(filepath + "m_init" + std::to_string(nz));
        // vti_writer_micro(state.Ms_field, mesh, filepath + "Ms");
        std::cout << "ncells= " << state.get_n_cells_() << std::endl;

        ////vti_writer_micro(state.m, mesh, filepath + "minit_nonnormalized" +
        /// std::to_string(nz));
        // state.m = normalize_handle_zero_vectors(state.m);
        // vti_writer_micro(state.m, mesh, filepath + "minit_renorm" +
        // std::to_string(nz));

        af::timer timer_llgterms = af::timer::start();
        af::array zee_field = af::constant(0, mesh.dims, f64);
        zee_field(af::span, af::span, af::span, 0) = zee;
        auto demag = LlgTerm(new DemagField(mesh));
        auto exch = LlgTerm(new ExchangeField(A));
        auto ext = LlgTerm(new ExternalField(zee_field));
        std::cout << "Llgterms assembled in " << af::timer::stop(timer_llgterms) << std::endl;

        LLGIntegrator llg(1, {demag, exch, ext});
        af::timer t_hys = af::timer::start();
        llg.relax(state, 1e-10);
        std::cout << "time minimize [af-s]: " << af::timer::stop(t_hys) << std::endl;
        // state.calc_mean_m_steps(stream,
        // afvalue(minimizer.llgterms_[minimizer.llgterms_.size() -
        // 1]->h(state)(0, 0, 0, 0))); vti_writer_micro(state.m, mesh, filepath
        // + "m_hysteresis_" + std::to_string(nz));
        state.write_vti(filepath + "m_relaxed" + std::to_string(nz));
        std::cout << nz << "\t" << state.meani(0) << "\t" << state.meani(1) << "\t" << state.meani(2) << "\t"
                  << afvalue(demag->h(state)(nx / 2, ny / 2, 0, 0)) * constants::mu0 << "\t" << std::endl;
        stream << nz << "\t" << state.meani(0) << "\t" << state.meani(1) << "\t" << state.meani(2) << "\t"
               << afvalue(demag->h(state)(nx / 2, ny / 2, 0, 0)) * constants::mu0 << "\t" << std::endl;
    }
    stream.close();
    return 0;
}
