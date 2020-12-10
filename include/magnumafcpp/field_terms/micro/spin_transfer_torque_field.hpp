#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"
#include "util/func.hpp"

namespace magnumafcpp {

class SpinTransferTorqueField : public MicroTerm {
  public:
    SpinTransferTorqueField(af::array polarization_field, double nu_dampinglike, double nu_field, double j_e);
    SpinTransferTorqueField(long int polarization_field_ptr, double nu_dampinglike, double nu_field, double j_e);

    WrappedArray polarization_field;

    double nu_dampinglike;
    double nu_fieldlike;
    double j_e;

  private:
    virtual af::array impl_H_in_Apm(const State& state) const override; ///< Effective field term contribution
};
} // namespace magnumafcpp
