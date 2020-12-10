#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

class ExchangeField : public MicroTerm {
  public:
    explicit ExchangeField(double A);
    explicit ExchangeField(af::array A_field);
    explicit ExchangeField(long int A_field_ptr);
    virtual af::array h(const State& state) const override;
    using MicroTerm::impl_E_in_J;
    virtual double impl_E_in_J(const State& state, const af::array& h) const override;

    double A{0}; //!< Exchange energy in [J/m]
    af::array A_field;

  private:
    // Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state) const;
};
} // namespace magnumafcpp
