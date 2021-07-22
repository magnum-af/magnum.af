#pragma once
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"
#include "state.hpp"

namespace magnumaf {

class AtomisticDipoleDipoleField : public AtomTerm {
  public:
    explicit AtomisticDipoleDipoleField(Mesh);

  private:
    virtual af::array impl_H_in_Apm(const State& state) const override;
    af::array Nfft;
};
} // namespace magnumaf
