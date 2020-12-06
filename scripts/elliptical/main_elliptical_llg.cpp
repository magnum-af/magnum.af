#include "arrayfire.h"
#include "magnum_af.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace magnumafcpp;

af::array zee_func(State state) {
    double field_Tesla = 0;
    double rate = 0.34e6;   //[T/s]
    double hzee_max = 0.25; //[T]
    if (state.t < hzee_max / rate)
        field_Tesla = rate * state.t;
    else if (state.t < 3 * hzee_max / rate)
        field_Tesla = -rate * state.t + 2 * hzee_max;
    else if (state.t < 4 * hzee_max / rate)
        field_Tesla = rate * state.t - 4 * hzee_max;
    else {
        field_Tesla = 0;
        std::cout << "WARNING ZEE time out of range" << std::endl;
    }
    array zee = constant(0.0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
    zee(span, span, span, 2) =
        constant(field_Tesla / state.constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
    return zee;
}

int main(int argc, char** argv) {
    std::string filepath(argc >= 1 ? argv[1] : "../Data/Testing");
    if (argc >= 1) {
        filepath.append("/");
    }
    if (argc >= 2) {
        setDevice(std::stoi(argv[2]));
    }
    std::string path_mrelax(argc > 3 ? argv[3] : "");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;
    std::cout.precision(24);
    info();

    // Parameter initialization
    Material material = Material();
    state.Ms = 2. / constants::mu0; //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    material.A = 1.5e-11;           //[J/m]
    material.Ku1 = 1.4e6;
    material.alpha = 0.02;

    const double x = 1000e-9, y = 6000e-9,
                 z = 5e-9; //[m] // Physical dimensions
    const int nx = 343;
    const int ny = 1920;
    const int nz = 2;

    // Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
    long int n_cells = 0; // Number of cells with Ms!=0
    State state(mesh, material, util::ellipse(mesh, n_cells));

    vti_writer_micro(state.Ms, mesh, (filepath + "Ms").c_str());
    vti_writer_micro(state.m, mesh, (filepath + "minit").c_str());
    std::cout << mesh << std::endl;

    std::vector<uptr_Fieldterm> llgterm;
    timer t_demag = af::timer::start();
    llgterm.push_back(uptr_Fieldterm(new DemagField(mesh, material)));
    std::cout << "Demag assembled in " << af::timer::stop(t_demag) << std::endl;
    llgterm.push_back(uptr_Fieldterm(new ExchangeField(mesh, material)));
    llgterm.push_back(uptr_Fieldterm(new UniaxialAnisotropyField(mesh, material)));
    LLGIntegrator Llg(llgterm);

    // Relaxation
    if (!exists(path_mrelax)) {
        Llg.relax(state, 1e-7);
        vti_writer_micro(state.m, mesh, (filepath + "mrelax").c_str());
        state.t = 0; // Setting t=0 for hysteresis
    } else {
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, mesh, path_mrelax);
    }

    std::cout << "n_cells= " << n_cells
              << ", should be a*b*M_PI*mesh.nz= " << mesh.nx / 2 * mesh.ny / 2 * M_PI * mesh.nz << std::endl;

    std::ofstream stream;
    stream.precision(12);
    stream.open((filepath + "m.dat").c_str());
    stream << "# t	<mx>" << std::endl;

    timer t_hys = af::timer::start();
    double rate = 0.34e6;                                          //[T/s]
    double hzee_max = 0.25;                                        //[T]
    Llg.llgterms.push_back(uptr_Fieldterm(new ExternalField(&zee_func))); // Rate in
                                                                   // T/s
    while (state.t < 4 * hzee_max / rate) {
        Llg.step(state);
        state.calc_mean_m(stream, n_cells, Llg.llgterms[Llg.llgterms.size() - 1]->h(state)(0, 0, 0, af::span));
        if (state.steps % 2000 == 0) {
            vti_writer_micro(state.m, mesh, (filepath + "m_hysteresis_" + std::to_string(state.steps)));
        }
    }

    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
