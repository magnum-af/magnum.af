#ifndef LBFGS_MINIMIZER_H
#define LBFGS_MINIMIZER_H
#include <memory>
#include <list>
#include <algorithm>
#include "arrayfire.h"
#include "../state.hpp"
#include "../llg_terms/LLGTerm.hpp"

//For second Method, use interface class:
//https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members

typedef std::shared_ptr<LLGTerm> LlgTerm; 
typedef std::vector<LlgTerm> LlgTerms; 

class LBFGS_Minimizer {
    public:
        LBFGS_Minimizer();

        void minimize(State&); // Minimization routine

        LlgTerms _llgterms;

        double get_time_h() const { return _time_h;};
    private:
        af::array _h(const State& m);// Effective Field 
        double _time_h{0};// Timer measuring calls to effective field _h
};

#endif
