#include "equations.hpp"
#include "field_terms/micro/demag_field.hpp"
#include "field_terms/micro/exchange_field.hpp"
#include "field_terms/micro/external_field.hpp"
#include "math.hpp"
#include "rk4.hpp"
#include "state.hpp"
#include "util/arg_parser.hpp"
#include "util/timer.hpp"
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <integrators/controller.hpp>
#include <integrators/rkf45.hpp>

using namespace magnumafcpp;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;

    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 100, ny = 25, nz = 1;
    const Mesh mesh{nx, ny, nz, x / nx, y / ny, z / nz};

    double alpha = 1.;
    const double Ms = 8e5;
    const double A = 1.3e-11;

    DemagField dmag(mesh, true, true, 0);
    ExchangeField exch(A);

    auto fieldterms = fieldterm::to_vec(dmag, exch);

    af::array m = af::constant(0, mesh::dims_v(mesh), f64);
    m(0, af::span, af::span, 1) = 1;
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(-1, af::span, af::span, 1) = 1;

    // Note, must capture alpha and filedterms by ref, as they are subject to change
    // Note state can be ref or not, as we overwrite changing values anyway
    State dummy_state(mesh, Ms, m);
    auto llg = [&alpha, &dummy_state, &fieldterms](const double t, const af::array& m_in) {
        dummy_state.t = t;
        dummy_state.m = m_in;
        const auto H_eff_in_Apm = fieldterm::Heff_in_Apm(fieldterms, dummy_state);
        return equations::LLG(alpha, dummy_state.m, H_eff_in_Apm);
    };

    std::ofstream os(outdir / "m.dat");
    os.precision(12);

    double t{0};
    double dt = 1.01e-15;
    Controller controller{};

    // relax
    StageTimer timer;
    while (t < 1e-9) {
        double error{0};
        af::array m_proposed;
        do {
            const auto [tn, yndy, rk_error] = RKF45(t, dt, m, llg);
            m_proposed = yndy;
            error = math::max_4d_abs(rk_error / controller.givescale(max(m, yndy)));
        } while (!controller.success(error, dt));
        t += dt;
        dt = controller.get_hnext();
        m = m_proposed;
        m = normalize(m);
        const auto [mx, my, mz] = math::mean_3d<double>(m);
        // std::cout << t << " " << mx << " " << my << " " << mz << std::endl;
        os << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
    }
    timer.print_stage("relax ");

    af::array external = af::constant(0.0, nx, ny, nz, 3, f64);
    external(af::span, af::span, af::span, 0) = -24.6e-3 / constants::mu0;
    external(af::span, af::span, af::span, 1) = +4.3e-3 / constants::mu0;
    fieldterms.push_back(fieldterm::to_uptr<ExternalField>(external));

    alpha = 0.02;
    while (t < 2e-9) {
        double error{0};
        af::array m_proposed;
        do {
            const auto [tn, yndy, rk_error] = RKF45(t, dt, m, llg);
            m_proposed = yndy;
            error = math::max_4d_abs(rk_error / controller.givescale(max(m, yndy)));
        } while (!controller.success(error, dt));
        t += dt;
        dt = controller.get_hnext();
        m = m_proposed;
        m = normalize(m);

        const auto [mx, my, mz] = math::mean_3d<double>(m);
        // std::cout << t << " " << mx << " " << my << " " << mz << std::endl;
        os << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
    }
    timer.print_stage("switch");
    timer.print_accumulated();
}
