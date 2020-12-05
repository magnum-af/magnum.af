#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"

namespace magnumafcpp {

class AtomisticDipoleDipoleField : public AtomisticTermBase {
  public:
    virtual af::array h(const State& state) const override;
    AtomisticDipoleDipoleField(Mesh);

  private:
    af::array Nfft;
};
} // namespace magnumafcpp
