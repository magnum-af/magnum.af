#pragma once
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"
#include "micro/external_field.hpp"

namespace magnumaf {

// TODO class AtomisticExternalField : public AtomTerm {
class AtomisticExternalField : public ExternalField {
  public:
    using ExternalField::ExternalField;

  private:
    virtual double impl_E_in_J(const State& state, const af::array& h) const override;
};

} // namespace magnumaf
