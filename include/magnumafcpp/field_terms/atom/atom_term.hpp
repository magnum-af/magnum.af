#pragma once
#include "arrayfire.h"
#include "field_term.hpp"
namespace magnumafcpp {

class AtomTerm : public FieldTerm {
  public:
    virtual ~AtomTerm() = default;

    ///< Calculating the atomistic energy Eex=-mu0/2 integral(M . Hex) dx
    using FieldTerm::E;
    virtual double E(const State& state, const af::array& h) const override {
        return -constants::mu0 / 2. * state.Ms *
               af::sum(af::sum(af::sum(af::sum(h * state.m, 0), 1), 2), 3).as(f64).scalar<double>();
    };
};

} // namespace magnumafcpp
