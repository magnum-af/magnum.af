#ifndef MINIMIZER_H
#define MINIMIZER_H
#include <memory>
#include "arrayfire.h"
#include "../state.hpp"
#include "../llg_terms/LLGTerm.hpp"


// Energy minimizer using a semi-implicit update scheme applying the Barzilian-Borwein (BB) rule for stepsize calculation.
typedef std::shared_ptr<LLGTerm> LlgTerm; 
typedef std::vector<LlgTerm> LlgTerms; 

class Minimizer {
    public:
        Minimizer(LlgTerms);
        Minimizer();

        LlgTerms llgterms;
    private:
      
};

#endif
