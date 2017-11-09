#ifndef LLGTerm_H
#define LLGTerm_H
#include "arrayfire.h"
#include "state.hpp"

class LLGTerm{
  public:
    virtual af::array h (const State& state) =0;
    virtual double E (const State& state)=0;
    virtual double get_cpu_time()=0;
};

#endif
