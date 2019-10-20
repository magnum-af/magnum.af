#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
#include "arrayfire.h"

namespace magnumaf{



class SpinTransferTorqueField : public LLGTerm {
  public:
    af::array h(const State& state); ///< Effective field term contribution
    float E(const State& state); ///< Energy contribution
    float E(const State& state, const af::array& h); ///< Micromagnetic energy for a given effective field

    float evaluation_timing{0.}; ///< Accumulated time used for evaluating all calls to h(const State&)
    float get_cpu_time(){return evaluation_timing;} ///< Get accumulated computation time for the all calls of h(const State&)

    SpinTransferTorqueField (af::array polarization_field, float nu_dampinglike, float nu_field, float j_e);
    SpinTransferTorqueField (long int polarization_field_ptr, float nu_dampinglike, float nu_field, float j_e);

    WrappedArray polarization_field;

    float nu_dampinglike;
    float nu_fieldlike;
    float j_e;
};
}// namespace magnumaf
