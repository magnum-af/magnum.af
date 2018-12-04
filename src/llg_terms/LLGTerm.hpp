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
    virtual double E (const State& state, const af::array& h)=0;///< Calculating the micromagnetic energy for a already calculated h field (to save computational cost)
    virtual double get_cpu_time()=0;
};

// Aliases used to initialize objects wich inherit from this class
using LlgTerm = std::shared_ptr<LLGTerm>; 
using LlgTerms = std::vector<LlgTerm>; 

#endif
