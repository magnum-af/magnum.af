#include "arg_parser.hpp"
#include "magnum_af.hpp"
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <iostream>

// Necessary template specialization for af::array
// Adapted form boost/numeric/odeint/algebra/vector_space_algebra.hpp
// and boost/numeric/odeint/external/vexcl/vexcl_norm_inf.hpp
namespace boost::numeric::odeint {
template <> struct vector_space_norm_inf<af::array> {
    typedef double result_type; // typedef definition is explicitly needed here for odeint internals.
    result_type operator()(const af::array& x) const {
        return af::max(af::max(af::max(af::max(af::abs(x), 0), 1), 2), 3).as(f64).scalar<double>();
    }
};
} // namespace boost::numeric::odeint

int main(int argc, char** argv) {
    std::cout << "Start" << std::endl;

    using namespace magnumafcpp;
    const auto [outdir_tmp, posargs] = ArgParser(argc, argv).outdir_posargs;
    const auto outdir = outdir_tmp;

    const double x_magaf = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 100, ny = 25, nz = 1;
    const Mesh mesh{nx, ny, nz, x_magaf / nx, y / ny, z / nz};

    double alpha = 1.;
    const double Ms = 8e5;
    const double A = 1.3e-11;

    auto dmag = DemagField(mesh, true, true, 0);
    auto exch = ExchangeField(A);
    auto fieldterms = fieldterm::mv_to_vec(dmag, exch);

    const auto type = af::dtype::f64;
    af::array m = af::constant(0, mesh::dims_v(mesh), type);
    m(0, af::span, af::span, 1) = 1;
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(-1, af::span, af::span, 1) = 1;

    State state(mesh, Ms, m);
    state.write_vti(outdir / "minit");

    auto llg = [&alpha, &state, &fieldterms](const af::array& m_in, af::array& dxdt, const double t) {
        // TODO: how often and where should we normalize?
        // This normalizes every dxdt evaluation, i.e. every function callback,  e.g. 7(!) times for DP54.
        // Meaning resulting m is not normalized, but input m is, leading to slight difference as in custom ode solvers.

        state.t = t;
        if (state.Ms_field.isempty()) {
            state.m = normalize(m_in);
        } else {
            state.m = normalize_handle_zero_vectors(m_in);
        }

        const auto H_eff_in_Apm = fieldterm::Heff_in_Apm(fieldterms, state);
        dxdt = equations::LLG(alpha, state.m, H_eff_in_Apm);
    };

    struct observe_m {
        std::filesystem::path outdir_;
        observe_m(std::filesystem::path outdir) : outdir_(outdir) {}

        void operator()(const af::array& m, double t) {
            const auto mean = af::mean(af::mean(af::mean(m, 0), 1), 2).as(f64);
            const auto mx = mean(0, 0, 0, 0).scalar<double>();
            const auto my = mean(0, 0, 0, 1).scalar<double>();
            const auto mz = mean(0, 0, 0, 2).scalar<double>();
            std::cout << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
            std::ofstream stream(outdir_ / "m.dat", std::ios::app);
            stream.precision(12);
            stream << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
        }
    };

    const double eps_abs = 1e-6;
    const double eps_rel = 1e-6;

    // Note:
    // constant: dt is first attempted step and output inverval combined, (could lead to slow startup when step is huge,
    // but is not a big deal, proably a timing issue)
    // times: dt is seperate from dt_view, i.e. first step dt not dependend on observer range
    // NOTE: Timings are not very reliable, const mode was once 1m40s, then 30s with fehlberg78!!!

    enum class intmode { constant, adaptive, times };
    // const auto mode = intmode::times;
    const auto mode = intmode::adaptive;
    // const auto mode = intmode::constant;

    std::ofstream stream(outdir / "m.dat");
    stream.precision(12);
    std::vector<std::size_t> steps;
    namespace ode = boost::numeric::odeint;

    // choosing an integrator via typedef
    // typedef ode::runge_kutta_fehlberg78<af::array, double, af::array, double, ode::vector_space_algebra>
    // stepper_type;
    typedef ode::runge_kutta_dopri5<af::array, double, af::array, double, ode::vector_space_algebra> stepper_type;
    // typedef ode::runge_kutta_cash_karp54<af::array, double, af::array, double, ode::vector_space_algebra>
    // stepper_type;
    // typedef ode::runge_kutta_fehlberg78<af::array, double, af::array, double, ode::vector_space_algebra>
    // stepper_type;
    const auto stepper = make_controlled(eps_abs, eps_rel, stepper_type());

    auto make_range = [](double start_time, double end_time, double dt) {
        std::vector<double> range;
        do {
            range.push_back(start_time);
            start_time += dt;
        } while (start_time <= end_time);
        return range;
    };

    auto integrate = [&](double start_time, double end_time, double dt, double dt_view) {
        if (mode == intmode::times) {
            const auto range1 = make_range(start_time, end_time, dt_view);
            // std::copy(std::begin(range1), std::end(range1), std::ostream_iterator<double>(std::cout, " "));
            steps.push_back(integrate_times(stepper, llg, m, range1, dt, observe_m{outdir}));
        } else if (mode == intmode::adaptive) {
            steps.push_back(integrate_adaptive(stepper, llg, m, start_time, end_time, dt, observe_m{outdir}));
        } else if (mode == intmode::constant) {
            steps.push_back(integrate_const(stepper, llg, m, start_time, end_time, dt_view, observe_m{outdir}));
        } else {
            std::cout << "ERROR" << std::endl;
        }
    };

    const double start_time = 0.0;
    const double middle_time = 1e-9;
    const double end_time = 2e-9;
    const double dt = 1e-11;                 // Initial step for ADS control
    const double dt_view = middle_time / 50; // Required time interval for dense output

    integrate(start_time, middle_time, dt, dt_view);

    std::cout << "steps=" << steps << std::endl;

    // Setting external field
    af::array external = af::constant(0.0, nx, ny, nz, 3, type);
    external(af::span, af::span, af::span, 0) = -24.6e-3 / constants::mu0;
    external(af::span, af::span, af::span, 1) = +4.3e-3 / constants::mu0;
    fieldterms.push_back(fieldterm::to_uptr<ExternalField>(external));
    alpha = 0.02;

    integrate(middle_time, end_time, dt, dt_view);

    std::cout << "steps=" << steps << std::endl;
    std::cout << "accum_steps=" << std::accumulate(steps.begin(), steps.end(), 0) << std::endl;
}
