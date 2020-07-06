#pragma once
#include "LLGTerm.hpp"
#include "arrayfire.h"
namespace magnumafcpp {

class AtomisticTermBase : public LLGTerm {
  public:
    virtual ~AtomisticTermBase(){};

    ///< Calculating the atomistic energy
    ///< Eex=-mu0/2 integral(M . Hex) dx
    double E(const State& state) {
        return -constants::mu0 / 2. * state.Ms *
               sum(sum(sum(sum(h(state) * state.m, 0), 1), 2), 3)
                   .scalar<double>();
    };
    ///< Calculating the atomistic energy for a already calculated h-field
    double E(const State& state, const af::array& h) {
        return -constants::mu0 / 2. * state.Ms *
               sum(sum(sum(sum(h * state.m, 0), 1), 2), 3).scalar<double>();
    };
};

} // namespace magnumafcpp
