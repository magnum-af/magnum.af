#include "integrators/llg_integrator.hpp"
#include "equations.hpp"
#include "math.hpp"
#include "state.hpp"
#include <fstream>
#include <memory>
#include <utility>

namespace magnumaf {

af::array LLGIntegrator::fheff(const State& state) const {
    af::timer timer_heff = af::timer::start();
    auto solution = fieldterm::Heff_in_Apm(llgterms, state);
    time_heff += af::timer::stop(timer_heff);
    return solution;
}

af::array LLGIntegrator::f(const State& state) const {
    // calls_fdmdt++;
    // timer_fdmdt=timer::start();
    if (dissipation_term_only) {
        if (alpha_field_) {
            return equations::LLG_damping(alpha_field_.value(), state.m, math::cross4(state.m, fheff(state)));
        } else {
            return equations::LLG_damping(alpha, state.m, math::cross4(state.m, fheff(state)));
        }

    } else {
        if (alpha_field_) {
            return equations::LLG(alpha_field_.value(), state.m, fheff(state));
        } else {
            return equations::LLG(alpha, state.m, fheff(state));
        }
    }
    // time_fdmdt+=af::timer::stop(timer_fdmdt);
}

// Energy calculation
double LLGIntegrator::E(const State& state) const { return fieldterm::Eeff_in_J(llgterms, state); }

void LLGIntegrator::relax(State& state, double precision, unsigned eval_E, unsigned iwritecout, bool verbose) {
    double start_time = state.t;
    af::timer t = af::timer::start();
    double E_prev = 1e20;
    while (std::fabs((E_prev - E(state)) / E_prev) > precision) {
        E_prev = E(state);
        for (unsigned i = 0; i < eval_E; i++) {
            step(state);
        }
        if (iwritecout > 0 and state.steps % iwritecout == 0) {
            if (verbose) {
                printf("relaxing: step=%llu, rel_diff= %e, <mx>=%f, <my>=%f, "
                       "<mz>=%f\n",
                       state.steps, std::fabs((E_prev - E(state)) / E_prev), state.meani(0), state.meani(1),
                       state.meani(2));
            }
        }
    }
    if (verbose) {
        printf("relaxed %e [s] with computation time of %e [af-s]. Current "
               "state.t= %e\n",
               state.t - start_time, af::timer::stop(t), state.t);
    }
}

void LLGIntegrator::integrate_dense(State& state, double time_in_s, double write_every_dt_in_s, std::ostream& os,
                                    bool verbose) {
    const double time_init = state.t;
    const double time_final = time_init + time_in_s;
    os << state << std::flush << '\n'; // flushing for pywrap

    double time_next_write = time_init + write_every_dt_in_s;
    while (state.t < time_final) {
        step(state);
        if (state.t >= time_next_write) {
            double dummy_err{0};
            double backstep = time_next_write - state.t; // val is negative
            if (time_next_write >= time_final - 1e-20) {
                // final backsep s.t. state.t = time_init + time_in_s;
                state.m += RKF45(state, backstep, dummy_err);
                state.t += backstep;
                os << state << std::flush << '\n';
                break; // otherwise while could perform an additional step due to rounding errors
            } else {
                State current_state = state;
                current_state.m += RKF45(current_state, backstep, dummy_err); // backwards step
                current_state.t += backstep;                                  // backwards
                os << current_state << std::flush << '\n';
                if (verbose) {
                    std::cout << "step: " << current_state << std::endl;
                }
            }
            time_next_write += write_every_dt_in_s;
            while (time_next_write < state.t) {
                std::cout << "LLGIntegrator:: Warning! AdaptiveRungeKutta already passed next time point, skipping "
                             "data point. Consider setting fewer data points. time_next_write < state.t: "
                          << time_next_write << " < " << state.t << std::endl;
                time_next_write += write_every_dt_in_s;
            }
        }
    }
}

void LLGIntegrator::integrate_dense(State& state, double time_in_s, double write_every_dt_in_s,
                                    const std::string& filename, bool verbose, bool append) {
    auto stream = std::ofstream(filename, append ? std::ios::app : std::ios::out);
    integrate_dense(state, time_in_s, write_every_dt_in_s, stream, verbose);
}

long int LLGIntegrator::h_addr(const State& state) const { return util::pywrap::send_copy_to_py(fheff(state)); }
} // namespace magnumaf
