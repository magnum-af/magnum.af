#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    using namespace af;
    std::cout << "argc" << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc > 1 ? argv[1] : "../Data/Testing");
    if (argc > 1)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    const double hzee_max = (argc > 3 ? std::stod(argv[3]) : 0.12); //[Tesla]
    const int quater_steps = (argc > 4 ? std::stoi(argv[4]) : 100);
    std::string path_mrelax(argc > 5 ? argv[3] : "");
    af::info();
    std::cout.precision(24);

    auto zee_func = [hzee_max, quater_steps](State state) -> af::array {
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
            constant(field_Tesla / constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
        return zee;
    };

    // Parameter initialization
    double Ms = 2. / constants::mu0; //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    double A = 1.5e-11;              //[J/m]
    double Ku1 = 1.4e6;

    const double x = 1000e-9, y = 6000e-9,
                 z = 5e-9; //[m] // Physical dimensions
    // const int nx = 343;
    // const int ny = 1920;
    // const int nz = 2;

    const int nx = 250;
    const int ny = 250;
    const int nz = 1;

    // Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);

    long int n_cells = 0; // Number of cells with Ms!=0
    State state(mesh, Ms, util::ellipse(mesh, n_cells));
    state.calc_mean_m(std::cout, n_cells);
    vti_writer_micro(state.m, mesh, (filepath + "minit_nonnormalized").c_str());
    vti_writer_micro(state.Ms_field, mesh, (filepath + "Ms").c_str());
    vti_writer_micro(state.m, mesh, (filepath + "minit").c_str());
    std::cout << mesh << std::endl;

    af::timer timer_llgterms = af::timer::start();
    // Minimizer minimizer("BB", 1e-10, 1e-5, 1e4, 100);
    LBFGS_Minimizer minimizer = LBFGS_Minimizer();
    minimizer.llgterms_.push_back(uptr_FieldTerm(new DemagField(mesh)));
    minimizer.llgterms_.push_back(uptr_FieldTerm(new ExchangeField(A)));
    minimizer.llgterms_.push_back(uptr_FieldTerm(new UniaxialAnisotropyField(Ku1)));
    std::cout << "llgterms assembled in " << af::timer::stop(timer_llgterms) << std::endl;

    // Relaxation
    if (!exists(path_mrelax)) {
        timer t = af::timer::start();
        minimizer.Minimize(state);
        std::cout << "timerelax [af-s]: " << af::timer::stop(t) << std::endl;
        vti_writer_micro(state.m, mesh, (filepath + "mrelax").c_str());
    } else {
        std::cout << "found mrelax. loading magnetization" << std::endl;
        vti_reader(state.m, state.mesh, path_mrelax);
    }

    std::ofstream stream;
    stream.precision(12);
    stream.open((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    state.calc_mean_m(stream, n_cells);

    timer t_hys = af::timer::start();
    double rate = hzee_max / quater_steps; //[T/s]
    minimizer.llgterms_.push_back(uptr_FieldTerm(new ExternalField(zee_func)));
    while (state.t < 4 * hzee_max / rate) {
        state.t += 1.;
        minimizer.Minimize(state);
        state.calc_mean_m(stream, afvalue(minimizer.llgterms_[minimizer.llgterms_.size() - 1]->H_in_Apm(state)(0, 0, 0, 2)));
        // TODO previously, maybe should be
        // considered//state.calc_mean_m(stream, n_cells,
        // afvalue(minimizer.llgterms_[minimizer.llgterms_.size()-1]->H_in_Apm(state)(0,
        // 0, 0, 2)));
        state.steps++;
        if (state.steps % 1 == 0) {
            vti_writer_micro(state.m, mesh, (filepath + "m_hysteresis_" + std::to_string(state.steps)).c_str());
        }
    }

    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
