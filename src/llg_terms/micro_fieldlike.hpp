#ifndef MICRO_FieldlikeTorque_H
#define MICRO_FieldlikeTorque_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../constants.hpp"
#include "../state.hpp"
#include "../func.hpp"
class FieldlikeTorque : public LLGTerm {
  public:
    af::array h(const State& state); ///< Effective field term contribution
    double E(const State& state); ///< Energy contribution
    double E(const State& state, const af::array& h); ///< Micromagnetic energy for a given effective field

    double evaluation_timing{0.}; ///< Accumulated time used for evaluating all calls to h(const State&)
    double get_cpu_time(){return evaluation_timing;} ///< Get accumulated computation time for the all calls of h(const State&)

    FieldlikeTorque (af::array polarization_field, double nu_field, double j_e);
    FieldlikeTorque (long int polarization_field_ptr, double nu_field, double j_e);

    af::array polarization_field;
    void set_polarization_field(long int aptr);
    long int get_polarization_field_addr();

    double nu_field;
    double j_e;
};
#endif
