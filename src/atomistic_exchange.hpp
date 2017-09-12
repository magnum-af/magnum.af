#ifndef ATOMISTIC_EXCHANGE_H
#define ATOMISTIC_EXCHANGE_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "state.hpp"
#include "func.hpp"

class ATOMISTIC_EXCHANGE : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double* get_cpu_time(){return &cpu_time;}

    ATOMISTIC_EXCHANGE (const Mesh& mesh);
    af::array filtr;

    double     cpu_time{0.};
    af::timer timer_solve;
    double     time_conv{0.};
    af::timer timer_conv;
    double     time_edges{0.};
    af::timer timer_edges;
};
#endif
