// Calculates pinned layer and reference layer hysteresis in x-axis
#include "field_terms/micro/demag_field.hpp"
#include "field_terms/micro/external_field.hpp"
#include "field_terms/micro/rkky_exchange_field.hpp"
#include "integrators/llg_integrator.hpp"
#include <cmath>

#include "util/arg_parser.hpp"

using namespace magnumaf;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;
    af::setBackend(AF_BACKEND_CPU);
    const double hzee_max = posargs.size() > 0 ? std::stod(posargs[0]) : 1.; //[Tesla]

    af::info();

    // Parameter initialization
    // z1: Pinned Layer
    // z2: Reference Layer
    const int nx = 1, ny = 1, nz = 2;
    const double dx = 2e-9;
    const unsigned steps_full_hysteresis = 100;

    // const double dx = 10e-9;
    auto _1D_field = af::dim4(nx, ny, nz, 1);

    const double RKKY_mJ_per_m2 = posargs.size() > 1 ? std::stod(posargs[1]) : -0.8;
    const double RKKY = RKKY_mJ_per_m2 * 1e-3 * dx;
    std::cout << "RKKY=" << RKKY << std::endl;

    const double Ms1 = posargs.size() > 2 ? std::stod(posargs[2]) : 1.4e6;
    std::cout << "Ms1=" << Ms1 << std::endl;
    const double Ms2 = 1e6;
    const double A = 15e-12; // Note: A is replaced by RKKY here

    // antoferomagnetic coupling simulated by external field
    // Jaf u1 . u2 A == Js H A dz
    // H = Jaf/(Js dz)
    // double Jaf = 0.72e-3;
    double Jaf = posargs.size() > 3 ? std::stod(posargs[3]) : 0.36e-3;
    std::cout << "Jaf=" << Jaf << std::endl;
    const double H_af = Jaf / (Ms1 * constants::mu0 * dx);
    std::cout << "H_af[Oe]=" << H_af << std::endl;
    std::cout << "H_af [T]=" << H_af * constants::mu0 << std::endl;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);

    // Initial magnetic field
    af::array m = af::constant(0.0, mesh::dims_v(mesh), f64);
    m(af::span, af::span, 0, 0) = 1.;
    m(af::span, af::span, 1, 0) = -1.;
    af::array Ms_field = af::constant(0.0, _1D_field, f64);
    Ms_field(af::span, af::span, 0) = Ms1;
    Ms_field(af::span, af::span, 1) = Ms2;
    State state(mesh, Ms_field, m);
    state.write_vti(outdir / "minit");

    auto rkky = uptr_FieldTerm(new RKKYExchangeField(RKKY_values(af::constant(RKKY, mesh::dims_v(mesh), f64)),
                                                     Exchange_values(af::constant(A, mesh::dims_v(mesh), f64)), mesh));

    auto demag = uptr_FieldTerm(new DemagField(mesh, true, true, 0));

    unsigned current_step = 0;
    // Defining H_zee via lamdas
    auto zee_func = [H_af, &current_step, hzee_max](const State& state) -> af::array {
        double field_Tesla = std::numeric_limits<double>::quiet_NaN();
        if (current_step < steps_full_hysteresis / 4) {
            field_Tesla = hzee_max * 4. * current_step / steps_full_hysteresis;
        } else if (current_step < 3 * steps_full_hysteresis / 4) {
            field_Tesla = -hzee_max * 4. * current_step / steps_full_hysteresis + 2 * hzee_max;
        } else if (current_step < steps_full_hysteresis) {
            field_Tesla = hzee_max * 4. * current_step / steps_full_hysteresis - 4 * hzee_max;
        } else {
            field_Tesla = 0;
            std::cout << "WARNING ZEE time out of range, setting external "
                         "field to zero"
                      << std::endl;
        }

        af::array zee = af::constant(0.0, mesh::dims_v(state.mesh), f64);
        zee(af::span, af::span, af::span, 0) = field_Tesla / constants::mu0;
        // Adding af-coupling field in pinned layer only:
        zee(af::span, af::span, 0, 0) += H_af; // simulated af coupling
        return zee;
    };

    auto external = uptr_FieldTerm(new ExternalField(zee_func));
    LLGIntegrator llg(1, {std::move(demag), std::move(rkky), std::move(external)});

    std::ofstream stream(outdir / "m.dat");
    stream.precision(12);

    std::vector<double> abs_my_rl; // Ref Layer my list

    for (unsigned i = 0; i <= steps_full_hysteresis; ++i) {
        current_step = i; // update for lambda
        llg.relax(state);
        const double Hx_component = llg.llgterms[2]->H_in_Apm(state)(0, 0, 1, 0).scalar<double>() * constants::mu0;
        const auto mx_z0 = state.m(0, 0, 0, 0).scalar<double>();
        const auto mx_z1 = state.m(0, 0, 1, 0).scalar<double>();
        const auto my_z0 = state.m(0, 0, 0, 1).scalar<double>();
        const auto my_z1 = state.m(0, 0, 1, 1).scalar<double>();
        abs_my_rl.push_back(std::abs(my_z1));

        std::cout << i << "\t" << Hx_component << "\t" << mx_z0 << "\t" << mx_z1 << "\t" << my_z0 << "\t" << my_z1
                  << std::endl;
        stream << i << "\t" << Hx_component << "\t" << mx_z0 << "\t" << mx_z1 << "\t" << my_z0 << "\t" << my_z1
               << std::endl;
        state.write_vti(outdir / ("m" + std::to_string(i)));
    }

    stream.close();
    stream.open(outdir / "plotfile.gpi");
    stream << "set terminal pdf;" << std::endl;
    stream << "set xlabel 'H_x [T]'" << std::endl;
    stream << "set ylabel 'm_x'" << std::endl;
    stream << "set output 'saf_hysteresis.pdf'" << std::endl;
    stream << "p 'm.dat' u 2:3 w l title 'm_x Pinned'";
    stream << ", '' u 2:4 w l title 'm_x Reference'";
    stream << ", '' u 2:6 w l title 'm_y Pinned'";
    stream << ", '' u 2:6 w l title 'm_y Reference'" << std::endl;
    // stream << "p 'm.dat' u 2:3 w l title 'Pinned'";
    // stream << ", '' u 2:4 w l title 'Reference'" << std::endl << std::endl;
    stream << "set terminal jpeg" << std::endl;
    stream << "set output 'saf_angle_error.jpg'" << std::endl;
    stream << "replot" << std::endl;
    stream.close();

    int syscall = std::system(("cd " + outdir.string() + " && gnuplot plotfile.gpi").c_str());
    if (syscall != 0) {
        std::cout << "syscall plotting with gnuplot failed" << std::endl;
    }
    return 0;
}
