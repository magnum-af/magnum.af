#include "equations.hpp"
#include "math.hpp"
#include "micro/demag_field.hpp"
#include "micro/exchange_field.hpp"
#include "micro/external_field.hpp"
#include "util/arg_parser.hpp"
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <utility>


namespace custom_ode {

template <class Value = double> struct rk54_fehlberg_coefficients_a1 : boost::array<Value, 1> {
    rk54_fehlberg_coefficients_a1() { (*this)[0] = static_cast<Value>(1) / static_cast<Value>(4); }
};

template <class Value = double> struct rk54_fehlberg_coefficients_a2 : boost::array<Value, 2> {
    rk54_fehlberg_coefficients_a2() {
        (*this)[0] = static_cast<Value>(3) / static_cast<Value>(32);
        (*this)[1] = static_cast<Value>(9) / static_cast<Value>(32);
    }
};

template <class Value = double> struct rk54_fehlberg_coefficients_a3 : boost::array<Value, 3> {
    rk54_fehlberg_coefficients_a3() {
        (*this)[0] = static_cast<Value>(1932) / static_cast<Value>(2197);
        (*this)[1] = static_cast<Value>(-7200) / static_cast<Value>(2197);
        (*this)[2] = static_cast<Value>(7296) / static_cast<Value>(2197);
    }
};

template <class Value = double> struct rk54_fehlberg_coefficients_a4 : boost::array<Value, 4> {
    rk54_fehlberg_coefficients_a4() {
        (*this)[0] = static_cast<Value>(439) / static_cast<Value>(216);
        (*this)[1] = static_cast<Value>(-8);
        (*this)[2] = static_cast<Value>(3680) / static_cast<Value>(513);
        (*this)[3] = static_cast<Value>(-845) / static_cast<Value>(4104);
    }
};

template <class Value = double> struct rk54_fehlberg_coefficients_a5 : boost::array<Value, 5> {
    rk54_fehlberg_coefficients_a5() {
        (*this)[0] = static_cast<Value>(-8) / static_cast<Value>(27);
        (*this)[1] = static_cast<Value>(2);
        (*this)[2] = static_cast<Value>(-3544) / static_cast<Value>(2565);
        (*this)[3] = static_cast<Value>(1859) / static_cast<Value>(4104);
        (*this)[4] = static_cast<Value>(-11) / static_cast<Value>(40);
    }
};

template <class Value = double> struct rk54_fehlberg_coefficients_b : boost::array<Value, 6> {
    rk54_fehlberg_coefficients_b() {
        (*this)[0] = static_cast<Value>(16) / static_cast<Value>(135);
        (*this)[1] = static_cast<Value>(0);
        (*this)[2] = static_cast<Value>(6656) / static_cast<Value>(12825);
        (*this)[3] = static_cast<Value>(28561) / static_cast<Value>(56430);
        (*this)[4] = static_cast<Value>(-9) / static_cast<Value>(50);
        (*this)[5] = static_cast<Value>(2) / static_cast<Value>(55);
    }
};

template <class Value = double> struct rk54_fehlberg_coefficients_db : boost::array<Value, 6> {
    rk54_fehlberg_coefficients_db() {
        (*this)[0] =
            static_cast<Value>(16) / static_cast<Value>(135) - static_cast<Value>(25) / static_cast<Value>(216);
        (*this)[1] = static_cast<Value>(0);
        (*this)[2] =
            static_cast<Value>(6656) / static_cast<Value>(12825) - static_cast<Value>(1408) / static_cast<Value>(2565);
        (*this)[3] =
            static_cast<Value>(28561) / static_cast<Value>(56430) - static_cast<Value>(2197) / static_cast<Value>(4104);
        (*this)[4] = static_cast<Value>(-9) / static_cast<Value>(50) - static_cast<Value>(-1) / static_cast<Value>(5);
        (*this)[5] = static_cast<Value>(2) / static_cast<Value>(55);
    }
};

template <class Value = double> struct rk54_fehlberg_coefficients_c : boost::array<Value, 6> {
    rk54_fehlberg_coefficients_c() {
        (*this)[0] = static_cast<Value>(0);
        (*this)[1] = static_cast<Value>(1) / static_cast<Value>(4);
        (*this)[2] = static_cast<Value>(3) / static_cast<Value>(8);
        (*this)[3] = static_cast<Value>(12) / static_cast<Value>(13);
        (*this)[4] = static_cast<Value>(1);
        (*this)[5] = static_cast<Value>(1) / static_cast<Value>(2);
    }
};

template <class State, class Value = double, class Deriv = State, class Time = Value,
          class Algebra = typename boost::numeric::odeint::algebra_dispatcher<State>::algebra_type,
          class Operations = typename boost::numeric::odeint::operations_dispatcher<State>::operations_type,
          class Resizer = boost::numeric::odeint::initially_resizer>
class runge_kutta_fehlberg54
    : public boost::numeric::odeint::explicit_error_generic_rk<6, 5, 5, 4, State, Value, Deriv, Time, Algebra,
                                                               Operations, Resizer> {

  public:
    typedef boost::numeric::odeint::explicit_error_generic_rk<6, 5, 5, 4, State, Value, Deriv, Time, Algebra,
                                                              Operations, Resizer>
        stepper_base_type;
    using state_type = typename stepper_base_type::state_type;
    using value_type = typename stepper_base_type::value_type;
    using deriv_type = typename stepper_base_type::deriv_type;
    using time_type = typename stepper_base_type::time_type;
    using algebra_type = typename stepper_base_type::algebra_type;
    using operations_type = typename stepper_base_type::operations_type;
    using resizer_typ = typename stepper_base_type::resizer_type;

    using stepper_type = typename stepper_base_type::stepper_type;
    using wrapped_state_type = typename stepper_base_type::wrapped_state_type;
    using wrapped_deriv_type = typename stepper_base_type::wrapped_deriv_type;

    explicit runge_kutta_fehlberg54(const algebra_type& algebra = algebra_type())
        : stepper_base_type(
              boost::fusion::make_vector(rk54_fehlberg_coefficients_a1<Value>(), rk54_fehlberg_coefficients_a2<Value>(),
                                         rk54_fehlberg_coefficients_a3<Value>(), rk54_fehlberg_coefficients_a4<Value>(),
                                         rk54_fehlberg_coefficients_a5<Value>()),
              rk54_fehlberg_coefficients_b<Value>(), rk54_fehlberg_coefficients_db<Value>(),
              rk54_fehlberg_coefficients_c<Value>(), algebra) {}
};
} // namespace custom_ode

namespace boost::numeric::odeint {
// Specializations for runge_kutta_cash_karp54
template <class State, class Value, class Deriv, class Time, class Algebra, class Operations, class Resize>
struct get_controller<custom_ode::runge_kutta_fehlberg54<State, Value, Deriv, Time, Algebra, Operations, Resize>> {
    using stepper_type = custom_ode::runge_kutta_fehlberg54<State, Value, Deriv, Time, Algebra, Operations, Resize>;
    using type = controlled_runge_kutta<stepper_type>;
};

} // namespace boost::numeric::odeint

// Necessary template specialization for af::array
// Adapted form boost/numeric/odeint/algebra/vector_space_algebra.hpp
// and boost/numeric/odeint/external/vexcl/vexcl_norm_inf.hpp
namespace boost::numeric::odeint {
template <> struct vector_space_norm_inf<af::array> {
    using result_type = double; // typedef definition is explicitly needed here for odeint internals.
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

    // m-normalizing LLG:
    auto llg_norm = [&alpha, &state, &fieldterms](const af::array& m_in, af::array& dxdt, const double t) {
        // TODO: how often and where should we normalize?
        // This normalizes every dxdt evaluation, i.e. every function callback,  e.g. 7(!) times for DP54.
        // Meaning resulting m is not normalized, but input m is, leading to slight difference as in custom ode solvers.

        state.t = t;
        if (state.Ms_field.isempty()) {
            state.m = util::normalize(m_in);
        } else {
            state.m = util::normalize_handle_zero_vectors(m_in);
        }

        const auto H_eff_in_Apm = fieldterm::Heff_in_Apm(fieldterms, state);
        dxdt = equations::LLG(alpha, state.m, H_eff_in_Apm);
    };

    // not-m-normalizing LLG:
    auto llg_regular = [&alpha, &state, &fieldterms](const af::array& m_in, af::array& dxdt, const double t) {
        state.t = t;
        state.m = m_in;
        const auto H_eff_in_Apm = fieldterm::Heff_in_Apm(fieldterms, state);
        dxdt = equations::LLG(alpha, state.m, H_eff_in_Apm);
    };

    struct observe_m {
        std::filesystem::path outdir_;
        explicit observe_m(std::filesystem::path outdir) : outdir_(std::move(outdir)) {}

        void operator()(const af::array& m, double t) {
            // Possible renorm location, when observer is called after every step, i.e. adaptive mode:
            // TODO // m = normalize(m);
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

    enum class IntegationMode { constant, adaptive, times, steps, steps_dense };
    // constexpr auto mode = IntegationMode::times;
    constexpr auto mode = IntegationMode::adaptive;
    // constexpr auto mode = IntegationMode::constant;
    // constexpr auto mode = IntegationMode::steps; // Exact same results ad adaptive (if same llg is used)
    // constexpr auto mode = IntegationMode::steps_dense; // Not exactly same as constant
    // // .., probably because dopri45 used high order interpolation

    std::ofstream stream(outdir / "m.dat"); // clears old file
    stream.precision(12);
    std::vector<std::size_t> steps;
    namespace ode = boost::numeric::odeint;

    // choosing an integrator via typedef
    // typedef ode::runge_kutta_fehlberg78<af::array, double, af::array, double, ode::vector_space_algebra>
    // stepper_type;
    using stepper_type = ode::runge_kutta_dopri5<af::array, double, af::array, double, ode::vector_space_algebra>;
    // typedef custom_ode::runge_kutta_fehlberg54<af::array, double, af::array, double, ode::vector_space_algebra>
    //     stepper_type;
    // typedef ode::runge_kutta_cash_karp54<af::array, double, af::array, double, ode::vector_space_algebra>
    // stepper_type;
    // typedef ode::runge_kutta_fehlberg78<af::array, double, af::array, double, ode::vector_space_algebra>
    // stepper_type;
    auto stepper = make_controlled(eps_abs, eps_rel, stepper_type()); // Note: cannot be const, if try_step is used.

    auto make_range = [](double start_time, double end_time, double dt) {
        std::vector<double> range;
        do {
            range.push_back(start_time);
            start_time += dt;
        } while (start_time <= end_time);
        return range;
    };

    auto integrate = [&](const double start_time, const double end_time, const double dt, const double dt_view) {
        switch (mode) {
        case IntegationMode::times: {
            const auto range = make_range(start_time, end_time, dt_view);
            // std::copy(std::begin(range), std::end(range), std::ostream_iterator<double>(std::cout, " "));
            steps.push_back(integrate_times(stepper, llg_norm, m, range, dt, observe_m{outdir}));
            break;
        }
        case IntegationMode::adaptive:
            steps.push_back(integrate_adaptive(stepper, llg_norm, m, start_time, end_time, dt, observe_m{outdir}));
            break;
        case IntegationMode::constant:
            steps.push_back(integrate_const(stepper, llg_norm, m, start_time, end_time, dt_view, observe_m{outdir}));
            break;
        case IntegationMode::steps: {
            double t = start_time;
            double try_dt = dt;
            auto observer = observe_m{outdir};
            observer(m, t);
            std::size_t sucessful_steps = 0;
            while (t < end_time) {
                enum class LlgNormalize { on, off };
                constexpr auto llg_normalize = LlgNormalize::on; // yields exactly the same as adaptive mode
                // constexpr auto llg_normalize = LlgNormalize::off; // different form adaptive mode

                // Note: careful: ode::controlled_step_result would map ...::fail to true and ...::success to false!
                // So we would need inverse logic:
                // bool step_was_not_sucess = true;
                // while (step_was_not_sucess) {
                //    step_was_not_sucess = stepper.try_step(llg_regular, m, t, try_dt);

                ode::controlled_step_result step_result;
                do {
                    if (t + try_dt > end_time) { // Assure, we end at end_time
                        try_dt = end_time - t;
                    }
                    switch (llg_normalize) {
                    case LlgNormalize::off:
                        step_result = stepper.try_step(llg_regular, m, t, try_dt);
                        break;
                    case LlgNormalize::on:
                        step_result = stepper.try_step(llg_norm, m, t, try_dt);
                        break;
                    }
                    // std::cout << "step: " << sucessful_steps << '\t' << t << '\t' << try_dt << std::endl;
                } while (step_result == ode::controlled_step_result::fail);

                sucessful_steps++;
                if (llg_normalize == LlgNormalize::off) {
                    m = util::normalize(m); // NOTE: using non_normalizing_llg
                }
                observer(m, t);
            }
            steps.push_back(sucessful_steps);
            break;
        }
        case IntegationMode::steps_dense: {
            double t = start_time;
            double try_dt = dt;
            auto observer = observe_m{outdir};
            observer(m, t);
            std::size_t sucessful_steps = 0;
            double t_next_observer_call = start_time + dt_view;
            while (t < end_time) {
                enum class LlgNormalize { on, off };
                constexpr auto llg_normalize = LlgNormalize::on;
                // constexpr auto llg_normalize = LlgNormalize::off;
                ode::controlled_step_result step_result;
                do {
                    // Call observer only at intervals dt_view
                    // TODO handle edge cases when step is much larger than view_dt
                    if (t + try_dt > t_next_observer_call) {
                        try_dt = t_next_observer_call - t;
                    }
                    if (t + try_dt > end_time) { // Assure, we end at end_time
                        try_dt = end_time - t;
                    }
                    switch (llg_normalize) {
                    case LlgNormalize::off:
                        step_result = stepper.try_step(llg_regular, m, t, try_dt);
                        break;
                    case LlgNormalize::on:
                        step_result = stepper.try_step(llg_norm, m, t, try_dt);
                        break;
                    }
                    // std::cout << "step: " << sucessful_steps << '\t' << t << '\t' << try_dt << std::endl;
                } while (step_result == ode::controlled_step_result::fail);

                sucessful_steps++;
                if (llg_normalize == LlgNormalize::off) {
                    m = util::normalize(m); // NOTE: using non_normalizing_llg
                }
                if (std::abs(t - t_next_observer_call) < 1e-16) { //  or t == end_time
                    // if (t == t_next_observer_call) { //  or t == end_time
                    observer(m, t);
                    t_next_observer_call += dt_view;
                }
            }
            steps.push_back(sucessful_steps);
            break;
        }
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
