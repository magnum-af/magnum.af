#include "stochastic_llg.hpp"
#include "equations.hpp"
#include "field_terms/field_term.hpp"
#include "util/util.hpp"
#include <memory>

namespace magnumaf {

// Energy calculation
double Stochastic_LLG::E(const State& state) const { return fieldterm::Eeff_in_J(fieldterms, state); }

af::array Stochastic_LLG::fheff(const State& state) const { return fieldterm::Heff_in_Apm(fieldterms, state); }

af::array Stochastic_LLG::detfdmdt(const State& state) const {
    fdmdt_calls++;
    return equations::LLG(util::DoubleOrArray(this->alpha), state.m, fheff(state));
}

af::array Stochastic_LLG::stochfdmdt(const State& state, const af::array& h_th) const {
    stochfdmdt_calls++;
    return equations::LLG(util::DoubleOrArray(this->alpha), state.m, fheff(state) + h_th);
}

} // namespace magnumaf
