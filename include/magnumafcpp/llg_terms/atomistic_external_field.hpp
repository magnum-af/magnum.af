#pragma once
#include "arrayfire.h"
#include "zee.hpp"

namespace magnumafcpp {

class AtomisticExternalField : public ExternalField {
  public:
    using ExternalField::ExternalField;
    double E(const State& state) override; // Energy contribution
    double
    E(const State& state,
      const af::array& h) override; ///< Calculating the micromagnetic energy
                                    ///< for a already calculated h field

  private:
};

} // namespace magnumafcpp
