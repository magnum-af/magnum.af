#pragma once
#include "../state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"

namespace magnumafcpp {

class AtomisticExchangeField : public AtomisticTermBase {
  public:
    AtomisticExchangeField(double J_atom);
    // Field contribution
    af::array h(const State& state);
    // CPU time
    double get_cpu_time() { return cpu_time; }

    const double J_atom; //!< Atomistic exchange energy [J]

    double cpu_time{0.};
    af::timer timer_solve;
    double time_conv{0.};
    af::timer timer_conv;
    double time_edges{0.};
    af::timer timer_edges;
};
} // namespace magnumafcpp
