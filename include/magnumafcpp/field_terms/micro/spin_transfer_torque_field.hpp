#pragma once
#include "util/func.hpp"
#include "state.hpp"
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"

namespace magnumafcpp {

class SpinTransferTorqueField : public IntegratorTermMeshBase {
  public:
    virtual af::array h(const State& state) const override; ///< Effective field term contribution

    SpinTransferTorqueField(af::array polarization_field, double nu_dampinglike, double nu_field, double j_e);
    SpinTransferTorqueField(long int polarization_field_ptr, double nu_dampinglike, double nu_field, double j_e);

    WrappedArray polarization_field;

    double nu_dampinglike;
    double nu_fieldlike;
    double j_e;
};
} // namespace magnumafcpp
