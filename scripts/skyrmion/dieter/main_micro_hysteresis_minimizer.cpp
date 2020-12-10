#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

void calc_mean_m(const State& state, std::ostream& myfile, double hzee) {
    array sum_dim3 = sum(sum(sum(state.m, 0), 1), 2);
    myfile << std::setw(12) << state.t << "\t" << afvalue(sum_dim3(span, span, span, 0)) << "\t"
           << afvalue(sum_dim3(span, span, span, 1)) << "\t" << afvalue(sum_dim3(span, span, span, 2)) << "\t" << hzee
           << std::endl;
}

double hzee_max = 2.;   //[T]
int quater_steps = 100; // One 4th of total steps

af::array zee_func(State state) {
    double field_Tesla = 0;
    double rate = hzee_max / quater_steps; //[T/s]
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
    zee(span, span, span, 0) =
        constant(field_Tesla / state.constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
    return zee;
}

int main(int argc, char** argv) {

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc > 1 ? argv[1] : "../Data/skyrmion_stoch");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();

    // Parameter initialization
    const int nx = 90, nz = 1;
    const double dx = 1.0e-9;
    const double dz = 0.6e-9;

    // Generating Objects
    Mesh mesh(nx, nx, nz, dx, dx, dz);
    Material material = Material();
    state.Ms = 580000;
    material.A = 15e-12;
    material.alpha = 1;
    material.D = 3e-3;
    material.Ku1 = 0.6e6;

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(span, span, span, 2) = -1;
    for (int ix = 0; ix < mesh.nx; ix++) {
        for (int iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > nx / 4.)
                m(ix, iy, span, 2) = 1.;
        }
    }

    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh, (filepath + "minit").c_str());

    // Relax
    af::timer timer_llgterms = af::timer::start();
    Minimizer minimizer("BB", 1e-10, 1e-5, 1e4, 100);
    // minimizer.llgterms.push_back( uptr_FieldTerm (new DemagField(mesh, material)));
    minimizer.llgterms.push_back(uptr_FieldTerm(new ExchangeField(mesh, material)));
    minimizer.llgterms.push_back(uptr_FieldTerm(new DmiField(mesh, material)));
    minimizer.llgterms.push_back(uptr_FieldTerm(new UniaxialAnisotropyField(mesh, material)));
    std::cout << "llgterms assembled in " << af::timer::stop(timer_llgterms) << std::endl;

    // obtaining relaxed magnetization
    timer t = af::timer::start();
    minimizer.minimize(state);
    std::cout << "timerelax [af-s]: " << af::timer::stop(t) << std::endl;
    vti_writer_micro(state.m, mesh, (filepath + "relax").c_str());

    // Hysteresis
    std::ofstream stream;
    stream.precision(12);
    stream.open((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    calc_mean_m(state, stream, 0);

    timer t_hys = af::timer::start();
    double rate = hzee_max / quater_steps; //[T/s]
    minimizer.llgterms.push_back(uptr_FieldTerm(new ExternalField(&zee_func)));
    while (state.t < 4 * hzee_max / rate) {
        minimizer.minimize(state);
        calc_mean_m(state, stream, afvalue(minimizer.llgterms[3]->H_in_Apm(state)(0, 0, 0, 0)));
        state.t += 1.;
        state.steps++;
        if (state.steps % 1 == 0) {
            vti_writer_micro(state.m, mesh, (filepath + "m_hysteresis_" + std::to_string(state.steps)).c_str());
        }
    }
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
