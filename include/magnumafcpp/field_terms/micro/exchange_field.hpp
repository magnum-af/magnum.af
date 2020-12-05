#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "field_terms/integrator_term_mesh_base.hpp"

namespace magnumafcpp {

class ExchangeField : public IntegratorTermMeshBase {
  public:
    ExchangeField(double A);
    ExchangeField(af::array A_field);
    ExchangeField(long int A_field_ptr);
    // Field contribution
    virtual af::array h(const State& state) const override;
    // Energy contribution
    // using LLGTerm::E;
    using IntegratorTermMeshBase::E; // bringing overload E(state, h) back into scope, hidden otherwise
    virtual double E(const State& state) const override;
    // TODO//virtual double E(const State& state, const af::array& h) const override;
    // IntegratorTermMeshBase::E(const State& state, const af::array& h) is called with h, not with h_withedges

    double A{0}; //!< Exchange energy in [J/m]
    af::array A_field;

  private:
    // Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state) const;
};
} // namespace magnumafcpp
