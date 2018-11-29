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

        double Minimize(State&);

        LlgTerms llgterms_;

        double GetTimeCalcHeff() const { return time_calc_heff_;}; ///< Accumulated time for calculation of Heff.
    private:
        af::array Gradient(const State&);///< Calculate gradient as energy-dissipation term of llg
        af::array Heff(const State& m);///< Effective Field 
        double E(const State&); ///< Calculate Energy
        double time_calc_heff_{0};///< Timer measuring calls to effective field _h
        int verbose_{3};///<TODO investigate definition, init value etc
        int maxIter_{10};///<TODO investigate definition, init value etc
        double mxmxhMax(const State& state);///< TODO investigate definition, init value etc
};

#endif
