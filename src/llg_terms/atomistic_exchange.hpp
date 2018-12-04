#ifndef ATOMISTIC_EXCHANGE_H
#define ATOMISTIC_EXCHANGE_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"

class ATOMISTIC_EXCHANGE : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    ATOMISTIC_EXCHANGE (const Mesh& mesh);

    double     cpu_time{0.};
    af::timer timer_solve;
    double     time_conv{0.};
    af::timer timer_conv;
    double     time_edges{0.};
    af::timer timer_edges;
};
#endif
