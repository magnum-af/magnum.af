#pragma once
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"
#include "state.hpp"

namespace magnumaf {

class AtomisticExchangeField : public AtomTerm {
  public:
    explicit AtomisticExchangeField(double J_atom);

  private:
    double J_atom; //!< Atomistic exchange energy [J]
    virtual af::array impl_H_in_Apm(const State& state) const override;
};
} // namespace magnumaf
