#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"

namespace magnumafcpp {

class AtomisticDipoleDipoleField : public AtomisticTermBase {
  public:
    af::array h(const State& state) const override;
    AtomisticDipoleDipoleField(Mesh);

  private:
    af::array Nfft;
};
} // namespace magnumafcpp
