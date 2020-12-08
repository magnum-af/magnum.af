#pragma once
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"
#include "micro/external_field.hpp"

namespace magnumafcpp {

// TODO class AtomisticExternalField : public AtomTerm {
class AtomisticExternalField : public ExternalField {
  public:
    using ExternalField::ExternalField;
    virtual double E(const State& state, const af::array& h) const override;

  private:
};

} // namespace magnumafcpp
