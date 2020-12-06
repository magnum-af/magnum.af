#pragma once
#include "arrayfire.h"
#include "atomistic_term_base.hpp"
#include "micro/external_field.hpp"

namespace magnumafcpp {

// TODO class AtomisticExternalField : public AtomisticTermBase {
class AtomisticExternalField : public ExternalField {
  public:
    using ExternalField::ExternalField;
    using LLGTerm::E; // Could be AtomisticTermBase::E if inherited from there
    virtual double E(const State& state, const af::array& h) const override;

  private:
};

} // namespace magnumafcpp
