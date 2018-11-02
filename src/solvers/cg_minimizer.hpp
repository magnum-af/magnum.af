#ifndef CG_MINIMIZER_H
#define CG_MINIMIZER_H
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

class CG_Minimizer {
    public:
        CG_Minimizer();

        af::array h(const State& m);// Effective Field 
        void minimize(State&); // Minimization routine

        LlgTerms llgterms;

        double get_time_h() const { return time_h;};
    private:
        double time_h{0};
        
};

#endif
