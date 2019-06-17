#ifndef MINIMIZER_H
#define MINIMIZER_H
#include <memory>
#include <list>
#include <algorithm>
#include "arrayfire.h"
#include "../func.hpp"
#include "../state.hpp"
#include "../llg_terms/LLGTerm.hpp"

//For second Method, use interface class:
//https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members

// Energy minimizer using a semi-implicit update scheme applying the Barzilian-Borwein (BB) rule for stepsize calculation.
typedef std::shared_ptr<LLGTerm> LlgTerm;
typedef std::vector<LlgTerm> LlgTerms;

class Minimizer {
    public:
        Minimizer(std::string scheme = "BB", double tau_min = 1e-10, double tau_max = 1e-5, double dm_max = 1e4, int samples = 10, bool info = false);

        af::array h(const State& m);// Effective Field
        void minimize(State&); // Minimization routine

        LlgTerms llgterms;

        double get_time_h() const { return time_h;};
    private:
        double time_h{0};
        af::array dm(const State& state);
        af::array m_next(const State& state, const double tau);

        const std::string scheme;
        const double tau_min;
        const double tau_max;
        const double dm_max;
        const unsigned int samples;
        const bool info; // if ture, prints step info to cout
        double E(const State& state);// only for testing, remove?

};

#endif
