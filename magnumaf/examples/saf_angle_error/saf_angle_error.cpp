// Calculates pinned layer and reference layer for a full 360 [deg] Hext-rotation
#include "field_terms/micro/demag_field.hpp"
#include "field_terms/micro/exchange_field.hpp"
#include "field_terms/micro/external_field.hpp"
#include "field_terms/micro/rkky_exchange_field.hpp"
#include "integrators/llg_integrator.hpp"
#include "util/arg_parser.hpp"

using namespace magnumaf;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;
    const double hzee_max = posargs.size() > 0 ? std::stod(posargs[0]) : 0.100; //[Tesla]
    // Parameter initialization
    // z1: Pinned Layer
    // z2: Reference Layer
    const int nx = 1, ny = 1, nz = 2;
    const double dx = 2e-9;

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

    //// Anisotropy
    //// const double K1 = 0.2e6;
    //// const double K2 = 0.2e6;
    // const double K1 = -1 / 2. * std::pow(Ms1, 2) * constants::mu0;
    // const double K2 = -1 / 2. * std::pow(Ms2, 2) * constants::mu0;
    // std::cout << "K1=" << K1 << std::endl;
    // std::cout << "K2=" << K2 << std::endl;
    // const std::array<double, 3> Ku1_axis = {0, 0, 1}; // TODO
    // af::array Ku1_field = af::constant(0.0, _1D_field, f64);
    // Ku1_field(af::span, af::span, 0) = K1;
    // Ku1_field(af::span, af::span, 1) = K2;
    // auto aniso = uptr_FieldTerm(new UniaxialAnisotropyField(Ku1_field, Ku1_axis));

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
        const double hx = hzee_max / constants::mu0 * std::cos(current_step * M_PI / 180.);
        const double hy = hzee_max / constants::mu0 * std::sin(current_step * M_PI / 180.);

        af::array zee = af::constant(0.0, mesh::dims_v(state.mesh), f64);
        zee(af::span, af::span, af::span, 0) = hx;
        zee(af::span, af::span, af::span, 1) = hy;
        zee(af::span, af::span, 0, 0) += H_af; // simulated af coupling
        return zee;
    };
    auto external = uptr_FieldTerm(new ExternalField(zee_func));
    LLGIntegrator llg(1.0, {std::move(demag), std::move(rkky), std::move(external)});
    // LLGIntegrator llg(1.0, {demag, rkky, external, aniso});

    std::ofstream stream(outdir / "m.dat");
    stream.precision(12);

    std::vector<double> abs_my_rl; // Ref Layer my list

    // for (unsigned i = 0; i < 360; i++) {
    for (unsigned i = 0; i <= 360; i += 20) {
        current_step = i;
        llg.relax(state);
        const double Hx_component = llg.llgterms[2]->H_in_Apm(state)(0, 0, 1, 0).scalar<double>() * constants::mu0;
        const auto my_z0 = state.m(0, 0, 0, 1).scalar<double>();
        const auto my_z1 = state.m(0, 0, 1, 1).scalar<double>();
        abs_my_rl.push_back(std::abs(my_z1));

        std::cout << i << "\t" << Hx_component << "\t" << my_z0 << "\t" << my_z1 << std::endl;
        stream << i << "\t" << Hx_component << "\t" << my_z0 << "\t" << my_z1 << std::endl;
    }
    stream.close();

    // for (auto it : abs_my_rl) {
    //    std::cout << it << std::endl;
    //}
    double sum = std::accumulate(abs_my_rl.begin(), abs_my_rl.end(), 0.0);
    double mean = sum / abs_my_rl.size();
    double max = *std::max_element(abs_my_rl.begin(), abs_my_rl.end());
    std::cout << "sum=" << sum << std::endl;
    std::cout << "max=" << max << std::endl;

    std::cout << "mean=" << mean << std::endl;
    stream.open(outdir / "table.dat");
    stream << "# dx <<  Ms1[J/T/m3] << RKKY[mJ/m2] << max(abs(my)) << J_af "
              "[J/m2]  << Haf[T] << mean(abs(my))"
           << std::endl;
    stream << dx << "\t" << Ms1 << "\t" << RKKY_mJ_per_m2 << "\t" << max << "\t" << Jaf << "\t" << H_af * constants::mu0
           << "\t" << mean << std::endl;
    stream.close();

    stream.open(outdir / "plotfile.gpi");
    stream << "set terminal pdf;" << std::endl;
    stream << "set xlabel 'H_x [T]'" << std::endl;
    stream << "set ylabel 'm_y'" << std::endl;
    stream << "set output 'saf_angle_error.pdf'" << std::endl;
    stream << "p 'm.dat' u 2:3 w l title 'Pinned'";
    stream << ", '' u 2:4 w l title 'Reference'" << std::endl << std::endl;
    stream << "set terminal jpeg" << std::endl;
    stream << "set output 'saf_angle_error.jpg'" << std::endl;
    stream << "replot" << std::endl;
    stream.close();

    int syscall = std::system(("cd " / outdir / " && gnuplot plotfile.gpi").c_str());
    if (syscall != 0) {
        std::cout << "syscall plotting with gnuplot failed" << std::endl;
    }

    return 0;
}
