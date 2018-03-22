#ifndef ZEE_H
#define ZEE_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "state.hpp"
#include "func.hpp"

class Zee : public LLGTerm {
  public:
    Zee(af::array zee_in);
    //TODO 
    Zee(double rate_in, double hzee_max_in);
    Zee(long int zee_in_addr);
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double get_cpu_time(){return cpu_time;}

    double rate;//[T/s]
    double hzee_max;//[T]
    af::array zee_field;
    double cpu_time{0.};
    af::timer timer;
};


#endif
