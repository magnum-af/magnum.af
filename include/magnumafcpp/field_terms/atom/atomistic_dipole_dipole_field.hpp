#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"

namespace magnumafcpp {

class AtomisticDipoleDipoleField : public AtomTerm {
  public:
    virtual af::array h(const State& state) const override;
    AtomisticDipoleDipoleField(Mesh);

  private:
    af::array Nfft;
};
} // namespace magnumafcpp
