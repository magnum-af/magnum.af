#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>
#include <fstream>

using namespace magnumafcpp;

int main(int argc, char **argv)
{
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
    const bool min_over_llg = false;

    const int nx = 250, ny = 250; //, nz = 1;
    const double x = 1000e-9, y = 1000e-9; //, z = 65e-9; //[m] // Physical dimensions
    const double Ms = 1.393e6; //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    const double A = 1.5e-11;  //[J/m]
    const double zee = 40e-3/constants::mu0;
    const double dz = 10e-9;//TODO set

    // Parameter initialization
    for (int nz = 2; nz < 22; nz++)
    {
        std::cout << "starting nz = " << nz << std::endl;

        //Generating Objects
        Mesh mesh(nx, ny, nz, x/nx, y/ny, dz);

        af::array m0 = mesh.init_vortex();
        m0(af::span, af::span, 0, af::span) = 0;
        State state(mesh, Ms, m0);
        state.write_vti((filepath + "m_init" + std::to_string(nz)).c_str());
        //vti_writer_micro(state.Ms_field, mesh, filepath + "Ms");
        std::cout << "ncells= " << state.get_n_cells_() << std::endl;

        af::timer timer_llgterms = af::timer::start();
        af::array zee_field = af::constant(0, mesh.dims, f64);
        zee_field(af::span, af::span, af::span, 0) = zee;
        auto demag = LlgTerm(new DemagField(mesh));
        //auto exch = LlgTerm(new SparseExchangeField(A, mesh));
        auto exch = LlgTerm(new ExchangeField(A));//only use with opencl !!!
        auto ext  = LlgTerm(new ExternalField(zee_field));
        std::cout << "Llgterms assembled in " << af::timer::stop(timer_llgterms) << std::endl;
        LBFGS_Minimizer minimizer = LBFGS_Minimizer({demag, exch, ext}, 1e-6, 1000, 0);
        minimizer.of_convergence.open(filepath + "minimizer_convergence.dat" + std::to_string(nz));
        LLGIntegrator llg(1, {demag, exch, ext});

        af::timer t_hys = af::timer::start();
        if (min_over_llg)
        {
            minimizer.Minimize(state);
        }
        else
        {
            llg.relax(state, 1e-10);
        }
        std::cout << "time minimize [af-s]: " << af::timer::stop(t_hys) << std::endl;
        state.write_vti((filepath + "m_relaxed" + std::to_string(nz)).c_str());
        af::array demagfield = demag->h(state);
        vti_writer_micro(demagfield, mesh, filepath + "relaxed_demag_nz" + std::to_string(nz));
        double demag_in_center = afvalue(demagfield(nx/2, ny/2, 0, 0)) * constants::mu0;
        double mean_demag = afvalue(af::mean( af::mean(demagfield(af::span, af::span, 0, 0), 0), 1)) * constants::mu0;
        std::cout << nz << "\t" << state.meani(0) << "\t" << state.meani(1) << "\t" << state.meani(2) << "\t" << demag_in_center << "\t" << mean_demag << "\t" << std::endl;
        stream    << nz << "\t" << state.meani(0) << "\t" << state.meani(1) << "\t" << state.meani(2) << "\t" << demag_in_center << "\t" << mean_demag << "\t" << std::endl;
    }
    stream.close();
    return 0;
}
