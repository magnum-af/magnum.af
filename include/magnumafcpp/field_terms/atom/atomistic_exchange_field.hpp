#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "field_terms/atom/atomistic_term.hpp"

namespace magnumafcpp {

class AtomisticExchangeField : public AtomisticTermBase {
  public:
    AtomisticExchangeField(double J_atom);
    virtual af::array h(const State& state) const override;

  private:
    double J_atom; //!< Atomistic exchange energy [J]
};
} // namespace magnumafcpp
