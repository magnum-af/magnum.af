#pragma once
#include "../state.hpp"
#include "arrayfire.h"
#include "integrator_term_mesh_base.hpp"
#include <array>

namespace magnumafcpp {

class CubicAnisotropyField : public IntegratorTermMeshBase {
  public:
    CubicAnisotropyField(double Kc1, double Kc2 = 0, double Kc3 = 0, std::array<double, 3> c1 = {1, 0, 0},
                         std::array<double, 3> c2 = {0, 1, 0});
    const double Kc1, Kc2, Kc3; // First, second and third order cubic anisotropy constants in [J/m^3]

    const std::array<double, 3> c1, c2, c3; // Pairwise orthogonal unit vectors.

    af::array h(const State& state);
    double get_cpu_time() { return 0; } // TODO//!< accumulated heff computation time in [s]
    double E(const State& state);
    double E(const State& state, const af::array& h);

  private:
    std::array<af::array, 3> h_1to3(const State& state);
};
} // namespace magnumafcpp
