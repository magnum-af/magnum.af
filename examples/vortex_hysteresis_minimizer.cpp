#include "field_terms/micro/demag_field.hpp"
#include "field_terms/micro/exchange_field.hpp"
#include "field_terms/micro/external_field.hpp"
#include "solvers/lbfgs_minimizer.hpp"
#include "util/arg_parser.hpp"
#include "util/geometry.hpp"
#include "util/timer.hpp"
#include "util/vtk_IO.hpp"
#include <filesystem>

namespace fs = std::filesystem;
void write_plotfile(fs::path outdir, std::string plotfile = "plotfile.gpi");
void plot(fs::path outdir, std::string plotfile = "plotfile.gpi");

using namespace magnumafcpp;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;
    write_plotfile(outdir);
    const double hzee_max = (posargs.size() > 0 ? std::stod(posargs[0]) : 0.12); //[Tesla]
    const unsigned int steps_full_hysteresis = (posargs.size() > 1 ? std::stoi(posargs[1]) : 200);

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
    vti_writer_micro(state.Ms_field, mesh, outdir / "Ms");
    std::cout << "ncells= " << state.get_n_cells_() << std::endl;

    vti_writer_micro(state.m, mesh, (outdir / "minit_nonnormalized").c_str());
    state.m = normalize_handle_zero_vectors(state.m);
    vti_writer_micro(state.m, mesh, (outdir / "minit_renorm").c_str());

    DemagField dmag(mesh);
    ExchangeField exch(A);
    ExternalField extr(zee_func);
    LBFGS_Minimizer minimizer =
        LBFGS_Minimizer(fieldterm::to_vec(std::move(dmag), std::move(exch), extr), 1e-6, 1000, 0);
    minimizer.of_convergence_.open(outdir / "minimizer_convergence.dat");

    std::ofstream stream;
    stream.precision(24);
    stream.open((outdir / "m.dat").c_str());
    const std::string oinfo("#<mx>    <my>    <mz>    Hz[T]");
    stream << oinfo << std::endl;
    std::cout << oinfo << std::endl;
    af::timer t_hys = af::timer::start();
    for (unsigned i = 0; i < steps_full_hysteresis; i++) {
        minimizer.Minimize(state);
        const auto extrHeff = extr.H_in_T(state).scalar<double>();
        const auto [mx, my, mz] = state.mean_m();
        std::cout << i << " " << mx << " " << my << " " << mz << ", Hx[T]=" << extrHeff << std::endl;
        stream << mx << " " << my << " " << mz << " " << extrHeff << std::endl;
        // stream << state << extrHeff << std::endl;
        if (state.steps % 10 == 0) {
            vti_writer_micro(state.m, mesh, (outdir / ("m_hysteresis_" + std::to_string(state.steps))).c_str());
        }
        state.steps++;
    }
    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    plot(outdir);
    return 0;
}

void write_plotfile(fs::path outdir, std::string plotfile) {
    std::ofstream o(outdir / plotfile);
    o << "set terminal pdf;\n"
      << "set xlabel 'mu_0*H_ext [T]';\n"
      << "set ylabel 'Average Magnetizaion';\n"
      << "set output 'hys.pdf';\n"
      << "p "
      << "'m.dat' u 5:2 w l t '<m_x>',"
      << "'m.dat' u 5:3 w l t '<m_y>',"
      << "'m.dat' u 5:4 w l t '<m_z>'";
}

void plot(fs::path outdir, std::string plotfile) {
    int syscall = std::system(("cd " + outdir.string() + " && gnuplot " + plotfile).c_str());
    syscall != 0 ? std::cout << "syscall plotting with gnuplot failed" << std::endl : std::cout << "";
}

void write_and_plot(fs::path outdir, std::string plotfile = "plotfile.gpi") {
    write_plotfile(outdir, plotfile);
    plot(outdir, plotfile);
}
