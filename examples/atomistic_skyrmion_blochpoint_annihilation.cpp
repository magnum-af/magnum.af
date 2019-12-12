#include "../src/magnum_af.hpp"
#include "arrayfire.h"

using namespace magnumafcpp;

typedef std::shared_ptr<LLGTerm> llgt_ptr;

int main(int argc, char **argv)
{

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    }

    std::string filepath(argc > 1 ? argv[1] : "../Data/skyrmion_stoch");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();

    // Parameter initialization
    const int nxy = 30, nz = 1;
    const double dx = 1e-9;

    //Generating Objects
    Mesh mesh(nxy, nxy, nz, dx, dx, dx);
    const double alpha = 1;
    // Micromagnetic parameters
    const double Ms = 1.1e6;
    const double A = 1.6e-11;
    const double D = 2 * 5.76e-3;
    const double Ku1 = 6.4e6;
    const double ext = 10e-3/constants::mu0;

    // Atomistic parameters
    const double J_atom = 2. * A * dx;
    const double D_atom = D * pow(dx, 2);
    const double K_atom = Ku1 * pow(dx, 3);
    const double p = Ms * pow(dx, 3); //Compensate nz=1 instead of nz=4

    // Initial magnetic field
    af::array m = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(af::span, af::span, af::span, 2) = -1;
    for (uint32_t ix = 0; ix < mesh.n0; ix++)
    {
        for (uint32_t iy = 0; iy < mesh.n1; iy++)
        {
            const double rx = double(ix) - mesh.n0 / 2.;
            const double ry = double(iy) - mesh.n1 / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > 30 / 4.)
                m(ix, iy, af::span, 2) = 1.;
        }
    }

    State state(mesh, p, m);
    vti_writer_atom(state.m, mesh, filepath + "minit");

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back(llgt_ptr(new AtomisticExchangeField(J_atom)));
    llgterm.push_back(llgt_ptr(new AtomisticDmiField(D_atom, {0, 0, -1})));
    llgterm.push_back(llgt_ptr(new AtomisticUniaxialAnisotropyField(K_atom, {0, 0, 1})));
    llgterm.push_back(llgt_ptr(new AtomisticDipoleDipoleField(mesh)));
    af::array zee = af::constant(0, mesh.dims, f64);
    zee(af::span, af::span, af::span, 2) = ext;
    llgterm.push_back(llgt_ptr(new AtomisticExternalField(zee)));

    LLGIntegrator Llg(alpha, llgterm);

    af::timer t = af::timer::start();
    Llg.relax(state, 1e-12, 100, 100);

    af::print("exch", af::mean(af::mean(af::mean(af::mean(llgterm[0]->h(state), 0), 1), 2), 3));
    af::print("dmi ", af::mean(af::mean(af::mean(af::mean(llgterm[1]->h(state), 0), 1), 2), 3));
    af::print("ani ", af::mean(af::mean(af::mean(af::mean(llgterm[2]->h(state), 0), 1), 2), 3));

    state.calc_mean_m(std::cout);
    vti_writer_atom(state.m, mesh, filepath + "relax");

    std::cout << "timerelax [af-s]: " << af::timer::stop(t) << ", steps = " << state.steps << std::endl;

    af::array last = constant(0, mesh.dims, f64);
    last(af::span, af::span, af::span, 2) = 1;

    std::vector<State> inputimages;
    inputimages.push_back(state);
    State state_last = state;
    state_last.m = last;
    inputimages.push_back(state_last);
    //inputimages.push_back(State(mesh, material, last));

    double n_interp = 60;
    double string_dt = 5e-14;
    const int string_steps = 10000;
    double string_abort_rel_diff = 1e-12;
    double string_abort_abs_diff = 1e-27;

    String string(state, inputimages, n_interp, string_dt, Llg);
    double barrier = string.run(filepath, string_abort_rel_diff, string_abort_abs_diff, string_steps);

    //without dmi and zee//double expected_barrier = 4.420526609492e-20;
    double expected_barrier = 4.420526609492e-20;//TODO
    double rel_diff_barrier = 2 * std::fabs(barrier - expected_barrier)/(barrier + expected_barrier);
    if(rel_diff_barrier < 1e-3)
    {
        std::cout << "Barrier as expected, relative difference to expected is " << rel_diff_barrier << std::endl;
    }
    else
    {
        std::cout << "Warning, barrier differs from expected with " << rel_diff_barrier << std::endl;
    }
    return 0;
}
