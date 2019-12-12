#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumafcpp{


class AtomisticExchangeField : public LLGTerm {
  public:
    AtomisticExchangeField(double J_atom);
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    AtomisticExchangeField ();
    const double J_atom{0.}; //!< Atomistic exchange energy [J]

    double     cpu_time{0.};
    af::timer timer_solve;
    double     time_conv{0.};
    af::timer timer_conv;
    double     time_edges{0.};
    af::timer timer_edges;
};
}// namespace magnumafcpp
