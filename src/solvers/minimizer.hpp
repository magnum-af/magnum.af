#pragma once
#include "../state.hpp"
#include "../llg_terms/LLGTerm.hpp"
#include "arrayfire.h"

namespace magnumaf{


//For second Method, use interface class:
//https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members

// Energy minimizer using a semi-implicit update scheme applying the Barzilian-Borwein (BB) rule for stepsize calculation.

class Minimizer {
    public:
        Minimizer(std::string scheme = "BB", float tau_min = 1e-10, float tau_max = 1e-5, float dm_max = 1e4, int samples = 10, bool info = false);

        af::array h(const State& m);// Effective Field
        void minimize(State&); // Minimization routine

        LlgTerms llgterms;

        float get_time_h() const { return time_h;};
    private:
        float time_h{0};
        af::array dm(const State& state);
        af::array m_next(const State& state, const float tau);

        const std::string scheme;
        const float tau_min;
        const float tau_max;
        const float dm_max;
        const unsigned int samples;
        const bool info; // if ture, prints step info to cout
        float E(const State& state);// only for testing, remove?

};

}// namespace magnumaf
