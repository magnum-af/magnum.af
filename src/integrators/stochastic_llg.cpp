#include "stochastic_llg.hpp"
#include "equations.hpp"
#include "field_terms/field_term.hpp"
#include "util/func.hpp"
#include <memory>

namespace magnumafcpp {

// Energy calculation
double Stochastic_LLG::E(const State& state) {
    double solution = 0.;
    for (unsigned i = 0; i < fieldterms.size(); ++i) {
        solution += fieldterms[i]->E(state);
    }
    return solution;
}

af::array Stochastic_LLG::fheff(const State& state) const {
    af::array solution = fieldterms[0]->h(state);
    for (unsigned i = 1; i < fieldterms.size(); ++i) {
        solution += fieldterms[i]->h(state);
    }
    return solution;
}

af::array Stochastic_LLG::detfdmdt(const State& state) const {
    fdmdt_calls++;
    return equations::LLG(alpha, state.m, fheff(state));
}

af::array Stochastic_LLG::stochfdmdt(const State& state, const af::array& h_th) const {
    stochfdmdt_calls++;
    return equations::LLG(alpha, state.m, fheff(state) + h_th);
}

} // namespace magnumafcpp
