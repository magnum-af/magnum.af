#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    std::cout << "argc" << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc > 1 ? argv[1] : "../Data/Testing");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath << std::endl;
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    const double hzee_max = (argc > 3 ? std::stod(argv[3]) : 0.12); //[Tesla]
    const unsigned int steps_full_hysteresis = (argc > 4 ? std::stoi(argv[4]) : 200);

    af::info();
    std::cout.precision(24);

    // Defining H_zee via lamdas
    auto zee_func = [hzee_max, steps_full_hysteresis](State state) -> af::array {
        double field_Tesla;
        if (state.steps < steps_full_hysteresis / 4) {
            field_Tesla = hzee_max * 4. * state.steps / steps_full_hysteresis;
        } else if (state.steps < 3 * steps_full_hysteresis / 4) {
            field_Tesla = -hzee_max * 4. * state.steps / steps_full_hysteresis + 2 * hzee_max;
        } else if (state.steps < steps_full_hysteresis) {
            field_Tesla = hzee_max * 4. * state.steps / steps_full_hysteresis - 4 * hzee_max;
        } else {
            field_Tesla = 0;
            std::cout << "WARNING ZEE time out of range, setting external "
                         "field to zero"
                      << std::endl;
        }
        // std::cout << "fild= "<< field_Tesla << std::endl;
        af::array zee = af::constant(0.0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
        zee(af::span, af::span, af::span, 0) =
            af::constant(field_Tesla / constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
        return zee;
    };

    // Parameter initialization
    const int nx = 250, ny = 250, nz = 1;
    const double x = 1600e-9, y = 1600e-9,
                 z = 65e-9; //[m] // Physical dimensions

    // Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
    double Ms = 1.393e6; //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    double A = 1.5e-11;  //[J/m]

    State state(mesh, Ms, util::init_vortex(mesh));
    vti_writer_micro(state.Ms_field, mesh, filepath + "Ms");
    std::cout << "ncells= " << state.get_n_cells_() << std::endl;

    vti_writer_micro(state.m, mesh, (filepath + "minit_nonnormalized").c_str());
    state.m = normalize_handle_zero_vectors(state.m);
    vti_writer_micro(state.m, mesh, (filepath + "minit_renorm").c_str());

    DemagField dmag(mesh);
    ExchangeField exch(A);
    ExternalField extr(zee_func);
    LBFGS_Minimizer minimizer =
        LBFGS_Minimizer(fieldterm::to_vec(std::move(dmag), std::move(exch), extr), 1e-6, 1000, 0);
    minimizer.of_convergence.open(filepath + "minimizer_convergence.dat");

    std::ofstream stream;
    stream.precision(24);
    stream.open((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    af::timer t_hys = af::timer::start();
    for (unsigned i = 0; i < steps_full_hysteresis; i++) {
        minimizer.Minimize(state);
        const auto extrHeff = extr.H_eff(state).scalar<double>();
        const auto [mx, my, mz] = state.mean_m();
        std::cout << "Step " << i << ": " << mx << " " << my << " " << mz << ", Hx[T]=" << constants::mu0 * extrHeff
                  << std::endl;
        stream << i << " " << mx << " " << my << " " << mz << " " << constants::mu0 * extrHeff << std::endl;
        // stream << state << extrHeff << std::endl;
        if (state.steps % 10 == 0) {
            vti_writer_micro(state.m, mesh, (filepath + "m_hysteresis_" + std::to_string(state.steps)).c_str());
        }
        state.steps++;
    }
    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
