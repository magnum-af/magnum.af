#pragma once
#include "arrayfire.h"
#include "atomistic_term_base.hpp"
#include "micro/external_field.hpp"

namespace magnumafcpp {

// TODO class AtomisticExternalField : public AtomisticTermBase {
class AtomisticExternalField : public ExternalField {
  public:
    using ExternalField::ExternalField;
    virtual double E(const State& state) const override;
    virtual double E(const State& state, const af::array& h) const override;

  private:
};

} // namespace magnumafcpp
