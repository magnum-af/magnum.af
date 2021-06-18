#include "equations.hpp"
#include "field_terms/micro/demag_field.hpp"
#include "field_terms/micro/exchange_field.hpp"
#include "field_terms/micro/external_field.hpp"
#include "math.hpp"
#include "rk4.hpp"
#include "state.hpp"
#include "util/arg_parser.hpp"
#include "util/timer.hpp"
#include "util/util.hpp" // normalize
#include <array>
#include <cmath>
#include <fstream>
#include <utility>

#include <gtest/gtest.h>

using namespace magnumafcpp;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;
    const std::size_t imax = 1000;
    const double dt = 1e-9 / imax;

    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 100, ny = 25, nz = 1;
    const Mesh mesh{nx, ny, nz, x / nx, y / ny, z / nz};

    double alpha = 1.;
    const double Ms = 8e5;
    const double A = 1.3e-11;

    DemagField dmag(mesh, true, true, 0);
    ExchangeField exch(A);

    auto fieldterms = fieldterm::to_vec(dmag, exch);

    auto f = [&alpha, mesh, Ms, &fieldterms](double t, af::array m_in) {
        State state(mesh, Ms, std::move(m_in));
        state.t = t;
        auto H_eff_in_Apm = fieldterm::Heff_in_Apm(fieldterms, state);
        return equations::LLG(alpha, state.m, H_eff_in_Apm);
    };

    double t = 0;
    af::array m = af::constant(0, mesh::dims_v(mesh), f64);
    m(0, af::span, af::span, 1) = 1;
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(-1, af::span, af::span, 1) = 1;

    std::ofstream os(outdir / "m.dat");
    os.precision(12);

    // relax
    StageTimer timer;
    for (std::size_t i = 0; i < imax; i++) {
        std::tie(t, m) = RK4(t, dt, m, f);
        m = util::normalize_handle_zero_vectors(m);
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
    for (std::size_t i = 0; i < imax; i++) {
        std::tie(t, m) = RK4(t, dt, m, f);
        m = util::normalize_handle_zero_vectors(m);
        const auto [mx, my, mz] = math::mean_3d<double>(m);
        // std::cout << t << " " << mx << " " << my << " " << mz << std::endl;
        os << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
    }
    timer.print_stage("switch");
    timer.print_accumulated();
}
