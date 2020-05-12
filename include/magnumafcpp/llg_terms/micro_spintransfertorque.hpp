#pragma once
#include "../func.hpp"
#include "../state.hpp"
#include "LLGTerm.hpp"
#include "arrayfire.h"

namespace magnumafcpp {

class SpinTransferTorqueField : public LLGTerm {
  public:
    af::array h(const State& state); ///< Effective field term contribution
    double E(const State& state);    ///< Energy contribution
    double
    E(const State& state,
      const af::array& h); ///< Micromagnetic energy for a given effective field

    double evaluation_timing{0.}; ///< Accumulated time used for evaluating all
                                  ///< calls to h(const State&)
    double get_cpu_time() {
        return evaluation_timing;
    } ///< Get accumulated computation time for the all calls of h(const State&)

    SpinTransferTorqueField(af::array polarization_field, double nu_dampinglike,
                            double nu_field, double j_e);
    SpinTransferTorqueField(long int polarization_field_ptr,
                            double nu_dampinglike, double nu_field, double j_e);

    WrappedArray polarization_field;

    double nu_dampinglike;
    double nu_fieldlike;
    double j_e;
};
} // namespace magnumafcpp
