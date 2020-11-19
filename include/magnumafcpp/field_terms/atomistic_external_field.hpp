#pragma once
#include "arrayfire.h"
#include "atomistic_term_base.hpp"
#include "external_field.hpp"

namespace magnumafcpp {

// TODO class AtomisticExternalField : public AtomisticTermBase {
class AtomisticExternalField : public ExternalField {
  public:
    using ExternalField::ExternalField;
    double E(const State& state) override; // Energy contribution
    double E(const State& state,
             const af::array& h) override; ///< Calculating the micromagnetic energy
                                           ///< for a already calculated h field

  private:
};

} // namespace magnumafcpp
