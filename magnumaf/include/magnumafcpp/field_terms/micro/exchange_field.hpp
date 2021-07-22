#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumaf {

class ExchangeField : public MicroTerm {
  public:
    explicit ExchangeField(double A);
    explicit ExchangeField(af::array A_field);
    explicit ExchangeField(long int A_field_ptr);
    double A{0}; //!< Exchange energy in [J/m]
    af::array A_field;

  private:
    // Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state) const;

    virtual af::array impl_H_in_Apm(const State& state) const override;
    virtual double impl_E_in_J(const State& state, const af::array& h) const override;
};
} // namespace magnumaf
