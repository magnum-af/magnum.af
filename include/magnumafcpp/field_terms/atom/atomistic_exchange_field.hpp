#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"

namespace magnumafcpp {

class AtomisticExchangeField : public AtomisticTermBase {
  public:
    AtomisticExchangeField(double J_atom);
    af::array h(const State& state) const override;

  private:
    double J_atom; //!< Atomistic exchange energy [J]
};
} // namespace magnumafcpp
