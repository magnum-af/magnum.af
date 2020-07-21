#pragma once
#include "../state.hpp"
#include "arrayfire.h"
#include "integrator_term_mesh_base.hpp"

namespace magnumafcpp {

class DmiField : public IntegratorTermMeshBase {
  public:
    DmiField(double D, std::array<double, 3> D_axis = {0, 0, 1});
    DmiField(af::array D_constants, std::array<double, 3> D_axis = {0, 0, 1});
    DmiField(const double D, double D_axis_x, double D_axis_y,
             double D_axis_z); //! wrapping only
    DmiField(long int D_constants_ptr, double D_axis_x, double D_axis_y,
             double D_axis_z); //! wrapping only

    af::array h(const State& state);
    double get_cpu_time() { return cpu_time; }

  private:
    double cpu_time{0.};
    const double D{0};
    const af::array D_constants;
    const std::array<double, 3> D_axis;
    void apply_boundary_condition(af::array& hfield, const State& state);
};
} // namespace magnumafcpp
