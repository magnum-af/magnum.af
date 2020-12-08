#include "solvers/cg_minimizer.hpp"
#include "util/func.hpp"
#include "util/misc.hpp"

namespace magnumafcpp {

void abort_on_empty_size(const vec_uptr_FieldTerm& fieldterms) {
    if (fieldterms.size() == 0) {
        std::cout << bold_red("ERROR: LBFGS_Minimizer::Heff: Number of "
                              "_llgterms == 0. Please add at least one term to "
                              "LBFGS_Minimizer.fieldterms! Aborting...")
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Calculation of effective field
af::array CG_Minimizer::Heff(const State& state) const {
    abort_on_empty_size(fieldterms);
    const auto solution = fieldterm::accumulate_Heff_in_Apm(fieldterms, state);
    return solution;
}

std::pair<double, af::array> CG_Minimizer::EnergyAndGradient(const State& state) const {
    abort_on_empty_size(fieldterms);
    const auto [heff, energy] = fieldterm::accumulate_Heff_in_Apm_and_E(fieldterms, state);
    return {energy, 1. / (constants::mu0 * state.Ms) * cross4(state.m, cross4(state.m, heff))};
}

void CG_Minimizer::Minimize(State& state) const {
    af::timer timer = af::timer::start();
    af::print("h in minimize", af::mean(Heff(state))); // TODO

    std::cout << "CG_Minimizer: time = " << af::timer::stop(timer) << std::endl;
}
} // namespace magnumafcpp
