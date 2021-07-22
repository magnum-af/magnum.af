#pragma once
#include "arrayfire.h"
#include "field_term.hpp"
namespace magnumaf {

class AtomTerm : public FieldTerm {
  public:
    virtual ~AtomTerm() = default;

  private:
    ///< Calculating the atomistic energy Eex=-mu0/2 integral(M . Hex) dx
    virtual double impl_E_in_J(const State& state, const af::array& h) const override {
        return -constants::mu0 / 2. * state.Ms *
               af::sum(af::sum(af::sum(af::sum(h * state.m, 0), 1), 2), 3).as(f64).scalar<double>();
    };
};

} // namespace magnumaf
