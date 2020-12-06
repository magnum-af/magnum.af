#pragma once
#include "LLGTerm.hpp"
#include "arrayfire.h"
namespace magnumafcpp {

class AtomisticTermBase : public LLGTerm {
  public:
    virtual ~AtomisticTermBase() = default;

    ///< Calculating the atomistic energy
    ///< Eex=-mu0/2 integral(M . Hex) dx
    virtual double E(const State& state) const override {
        return -constants::mu0 / 2. * state.Ms *
               af::sum(af::sum(af::sum(af::sum(h(state) * state.m, 0), 1), 2), 3).as(f64).scalar<double>();
    };
    ///< Calculating the atomistic energy for a already calculated h-field
    virtual double E(const State& state, const af::array& h) const override {
        return -constants::mu0 / 2. * state.Ms *
               af::sum(af::sum(af::sum(af::sum(h * state.m, 0), 1), 2), 3).as(f64).scalar<double>();
    };
};

} // namespace magnumafcpp
