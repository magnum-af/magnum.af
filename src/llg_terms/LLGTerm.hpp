#ifndef LLGTerm_H
#define LLGTerm_H
#include <memory>
#include "arrayfire.h"
#include "../state.hpp"
// Abstract basis class for all terms in the LLG equation.
class LLGTerm{
  public:
    virtual af::array h (const State& state) =0;
    virtual double E (const State& state)=0;
    virtual double get_cpu_time()=0;
};

// Typedefs used to initialize objects wich inherit from this class
typedef std::shared_ptr<LLGTerm> LlgTerm; 
typedef std::vector<LlgTerm> LlgTerms; 

#endif
