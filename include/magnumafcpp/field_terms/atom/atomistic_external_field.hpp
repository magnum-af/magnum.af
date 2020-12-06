#pragma once
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"
#include "micro/external_field.hpp"

namespace magnumafcpp {

// TODO class AtomisticExternalField : public AtomTerm {
class AtomisticExternalField : public ExternalField {
  public:
    using ExternalField::ExternalField;
    using FieldTerm::E; // Could be AtomTerm::E if inherited from there
    virtual double E(const State& state, const af::array& h) const override;

  private:
};

} // namespace magnumafcpp
