#pragma once
#include "../state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"

namespace magnumafcpp {

class AtomisticDipoleDipoleField : public AtomisticTermBase {
  public:
    // Field contribution
    af::array h(const State& state);
    // CPU time
    double get_cpu_time() { return cpu_time; }

    AtomisticDipoleDipoleField(Mesh);

    af::array Nfft;
    double cpu_time{0.};
    af::timer timer_demagsolve;
};
} // namespace magnumafcpp
