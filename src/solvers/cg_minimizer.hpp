#pragma once
#include "arrayfire.h"
#include "../state.hpp"
#include "../llg_terms/LLGTerm.hpp"

namespace magnumaf{


//For second Method, use interface class:
//https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members

//! Preconditioned non-linear conjugate gradient minimizer.

//! Preconditioned Non-Linear Conjugate Gradient Minimizer (courtesy of Thomas Schrefl et.al.)
//!
class CG_Minimizer {
    public:
        CG_Minimizer();
        void Minimize(State&); // Minimization routine
        LlgTerms llgterms_;
        double GetTimeCalcHeff() const { return time_calc_heff_;}; ///< Accumulated time for calculation of Heff.
    private:
        af::array Heff(const State& m);///< Effective Field
        double EnergyAndGradient(const State& state, af::array& gradient);
        double time_calc_heff_{0};///< Timer measuring calls to effective field _h
};

}// namespace magnumaf
