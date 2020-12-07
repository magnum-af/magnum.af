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

void sync(const State& state) {
    if (state.afsync) {
        af::sync();
    }
}

// Calculation of effective field
af::array CG_Minimizer::Heff(const State& state) const {
    abort_on_empty_size(fieldterms);
    af::timer timer = af::timer::start();
    af::array solution = fieldterms[0]->h(state);
    for (unsigned i = 1; i < fieldterms.size(); ++i) {
        solution += fieldterms[i]->h(state);
    }
    sync(state);
    time_calc_heff_ += af::timer::stop(timer);
    return solution;
}

std::pair<double, af::array> CG_Minimizer::EnergyAndGradient(const State& state) const {
    abort_on_empty_size(fieldterms);
    af::timer timer = af::timer::start();
    // Avoiding array with zeros, starting loop with second term in llgterms
    af::array h = fieldterms[0]->h(state);
    double energy = fieldterms[0]->E(state, h);
    for (unsigned i = 1; i < fieldterms.size(); ++i) {
        af::array temp_h = fieldterms[i]->h(state);
        h += temp_h;
        energy += fieldterms[i]->E(state, temp_h);
    }
    sync(state);
    time_calc_heff_ += af::timer::stop(timer);
    return {energy, 1. / (constants::mu0 * state.Ms) * cross4(state.m, cross4(state.m, h))};
}

void CG_Minimizer::Minimize(State& state) const {
    af::timer timer = af::timer::start();
    af::print("h in minimize", af::mean(Heff(state))); // TODO

    std::cout << "CG_Minimizer: time = " << af::timer::stop(timer) << std::endl;
}
} // namespace magnumafcpp
