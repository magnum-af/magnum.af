#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

class DmiField : public MicroTerm {
  public:
    DmiField(double D, std::array<double, 3> D_axis = {0, 0, 1});
    DmiField(af::array D_constants, std::array<double, 3> D_axis = {0, 0, 1});
    DmiField(const double D, double D_axis_x, double D_axis_y,
             double D_axis_z); //! wrapping only
    DmiField(long int D_constants_ptr, double D_axis_x, double D_axis_y,
             double D_axis_z); //! wrapping only

  private:
    double D{0};
    af::array D_constants;
    std::array<double, 3> D_axis;
    void apply_boundary_condition(af::array& hfield, const State& state) const;
    virtual af::array impl_H_in_Apm(const State& state) const override;
};
} // namespace magnumafcpp
