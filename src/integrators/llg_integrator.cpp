#include "integrators/llg_integrator.hpp"
#include "equations.hpp"
#include "state.hpp"
#include "util/func.hpp"
#include <memory>

namespace magnumafcpp {

LLGIntegrator::LLGIntegrator(double alpha, std::string scheme, Controller controller, bool dissipation_term_only)
    : AdaptiveRungeKutta(scheme, controller), alpha(alpha), dissipation_term_only(dissipation_term_only) {}

LLGIntegrator::LLGIntegrator(double alpha, vec_uptr_FieldTerm llgterms, std::string scheme, Controller controller,
                             bool dissipation_term_only)
    : AdaptiveRungeKutta(scheme, controller), alpha(alpha), llgterms(std::move(llgterms)),
      dissipation_term_only(dissipation_term_only) {}

LLGIntegrator::LLGIntegrator(double alpha, std::initializer_list<movable_il<uptr_FieldTerm>> il, std::string scheme,
                             Controller controller, bool dissipation_term_only)
    : LLGIntegrator(alpha, {std::make_move_iterator(std::begin(il)), std::make_move_iterator(std::end(il))}, scheme,
                    controller, dissipation_term_only) {}

af::array LLGIntegrator::fheff(const State& state) const {
    af::timer timer_heff = af::timer::start();
    const auto solution = fieldterm::accumulate_Heff_in_Apm(llgterms, state);
    time_heff += af::timer::stop(timer_heff);
    return solution;
}

af::array LLGIntegrator::f(const State& state) const {
    // calls_fdmdt++;
    // timer_fdmdt=timer::start();
    if (dissipation_term_only) {
        return equations::LLG_damping(alpha, state.m, cross4(state.m, fheff(state)));

    } else {
        return equations::LLG(alpha, state.m, fheff(state));
    }
    // time_fdmdt+=af::timer::stop(timer_fdmdt);
}

// Energy calculation
double LLGIntegrator::E(const State& state) const { return fieldterm::accumulate_E_in_J(llgterms, state); }

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

long int LLGIntegrator::h_addr(const State& state) const {
    af::array* heff = new af::array(fheff(state));
    return (long int)heff->get();
}
} // namespace magnumafcpp
