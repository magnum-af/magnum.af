#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"
namespace magnumaf {

class DMI_D2d_Field : public MicroTerm {
  public:
    explicit DMI_D2d_Field(double D_in_J_per_m2, bool PBC = false) : D_in_J_per_m2_(D_in_J_per_m2), PBC_(PBC) {}
    double D_in_J_per_m2_{}; ///< DMI constant in [J/m2]
    bool PBC_{};             ///< if true, enables periodic boundary conditions (PBC)

  private:
    virtual af::array impl_H_in_Apm(const State& state) const override;
};

} // namespace magnumaf
