#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"

namespace magnumafcpp {

class ExchangeField : public IntegratorTermMeshBase {
  public:
    ExchangeField(double A);
    ExchangeField(af::array A_field);
    ExchangeField(long int A_field_ptr);
    virtual af::array h(const State& state) const override;
    using IntegratorTermMeshBase::E;
    virtual double E(const State& state, const af::array& h) const override;

    double A{0}; //!< Exchange energy in [J/m]
    af::array A_field;

  private:
    // Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state) const;
};
} // namespace magnumafcpp
