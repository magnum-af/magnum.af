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
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <integrators/controller.hpp>
#include <integrators/rkf45.hpp>

using namespace magnumaf;

template <typename T = double, typename Y = af::array> class RkIntegrator {
  public:
    explicit RkIntegrator(T dt = 1.01e-15) : dt_(dt) {}
    T t_{0};
    T dt_{1.01e-15};
    Controller controller{};
    // const bool normalize_ = true;

    template <typename F, typename... Args> void stepRKF45(Y& yn, F f, Args... args) {
        Y y_proposed_;
        T error{0};
        do {
            const auto [tn, yndy, rk_error] = RKF45(t_, dt_, yn, f, args...);
            y_proposed_ = yndy;
            error = math::max_4d_abs(rk_error / controller.givescale(max(yn, yndy)));
        } while (!controller.success(error, dt_));
        t_ += dt_;
        dt_ = controller.get_hnext();
        yn = util::normalize(y_proposed_);
        // handle at callsite// yn = util::normalize(m);
    }

    template <typename F, typename... Args>
    void integrateRKF45(const T time, Y& yn, F f, std::ostream& os, Args... args) {
        const T current_time = t_;
        while (t_ < current_time + time) {
            stepRKF45(yn, f, args...);
            const auto [mx, my, mz] = math::mean_3d<double>(yn);
            os << t_ << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
        }
    }

    template <typename F, typename... Args> void stepRK4(Y& yn, F f, Args... args) {
        std::tie(t_, yn) = RK4(t_, dt_, yn, f, args...);
        yn = util::normalize(yn);
    }
};

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

    // relax
    StageTimer timer;

    RkIntegrator rki{};
    // RK_Integrator rki{1e-12};

    rki.integrateRKF45(1e-9, m, llg, os);

    // while (rki.t_ < 1e-9) {
    //    // rki.stepRK4(m, llg);
    //    rki.stepRKF45(m, llg);
    //    const auto [mx, my, mz] = math::mean_3d<double>(m);
    //    // std::cout << t << " " << mx << " " << my << " " << mz << std::endl;
    //    os << rki.t_ << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
    //}

    timer.print_stage("relax ");

    af::array external = af::constant(0.0, nx, ny, nz, 3, f64);
    external(af::span, af::span, af::span, 0) = -24.6e-3 / constants::mu0;
    external(af::span, af::span, af::span, 1) = +4.3e-3 / constants::mu0;
    fieldterms.push_back(fieldterm::to_uptr<ExternalField>(external));

    alpha = 0.02;

    rki.integrateRKF45(0.997e-9, m, llg, os);

    // while (rki.t_ < 2e-9) {
    //    // rki.stepRK4(m, llg);
    //    rki.stepRKF45(m, llg);
    //    const auto [mx, my, mz] = math::mean_3d<double>(m);
    //    // std::cout << t << " " << mx << " " << my << " " << mz << std::endl;
    //    os << rki.t_ << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
    //}

    timer.print_stage("switch");
    timer.print_accumulated();
}
