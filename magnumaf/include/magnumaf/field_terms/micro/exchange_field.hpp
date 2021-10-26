#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumaf {

class ExchangeField : public MicroTerm {
  public:
    explicit ExchangeField(double A) : A(A) {}
    double A{}; //!< Exchange energy in [J/m]

  private:
    virtual af::array impl_H_in_Apm(const State& state) const override;
    virtual double impl_E_in_J(const State& state, const af::array& h) const override;
};
} // namespace magnumaf
