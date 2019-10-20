#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumaf{


class AtomisticExchangeField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    float E(const State& state);
    float E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    float get_cpu_time(){return cpu_time;}

    AtomisticExchangeField (const Mesh& mesh);

    float     cpu_time{0.};
    af::timer timer_solve;
    float     time_conv{0.};
    af::timer timer_conv;
    float     time_edges{0.};
    af::timer timer_edges;
};
}// namespace magnumaf
