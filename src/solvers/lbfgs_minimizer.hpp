#ifndef LBFGS_MINIMIZER_H
#define LBFGS_MINIMIZER_H
#include <memory>
#include <list>
#include <algorithm>
#include "arrayfire.h"
#include "../state.hpp"
#include "../func.hpp"
#include "../llg_terms/LLGTerm.hpp"

//For second Method, use interface class:
//https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members

class LBFGS_Minimizer {
    public:
        LBFGS_Minimizer();

        void Minimize(State&); // Minimization routine

        LlgTerms llgterms_;

        double GetTimeCalcHeff() const { return time_calc_heff_;};
    private:
        af::array Gradient(const State&);// Calculate gradient as energy-dissipation term of llg
        af::array Heff(const State& m);// Effective Field 
        double time_calc_heff_{0};// Timer measuring calls to effective field _h
};

#endif
