#ifndef MICRO_DampinglikeTorque_H
#define MICRO_DampinglikeTorque_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../constants.hpp"
#include "../state.hpp"
#include "../func.hpp"
class DampinglikeTorque : public LLGTerm {
  public:
    af::array h(const State& state); ///< Effective field term contribution
    double E(const State& state); ///< Energy contribution
    double E(const State& state, const af::array& h); ///< Micromagnetic energy for a given effective field

    double evaluation_timing{0.}; ///< Accumulated time used for evaluating all calls to h(const State&)
    double get_cpu_time(){return evaluation_timing;} ///< Get accumulated computation time for the all calls of h(const State&)

    DampinglikeTorque (af::array polarization_field, double nu_damp, double j_e);

    af::array polarization_field;
    double nu_damp;
    double j_e;
};
#endif
