#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "double_or_array.hpp"
#include "integrator_term_mesh_base.hpp"
#include "unit_vector_or_array.hpp"
#include <array>

namespace magnumafcpp {

class CubicAnisotropyField : public IntegratorTermMeshBase {
  public:
    CubicAnisotropyField(double Kc1, double Kc2 = 0, double Kc3 = 0, std::array<double, 3> c1 = {1, 0, 0},
                         std::array<double, 3> c2 = {0, 1, 0});
    CubicAnisotropyField(af::array Kc1_array, af::array Kc2_array, af::array Kc3_array, af::array c1_array,
                         af::array c2_array);
    // Wrapping only
    CubicAnisotropyField(double Kc1, double Kc2, double Kc3, double c1x, double c1y, double c1z, double c2x, double c2y,
                         double c2z);
    CubicAnisotropyField(long int Kc1_array_ptr, long int Kc2_array_ptr, long int Kc3_array_ptr, long int c1_array_ptr,
                         long int c2_array_ptr);

    const DoubleOrArray Kc1, Kc2, Kc3; // First, second and third order cubic anisotropy constants in [J/m^3]

    /// Pairwise orthogonal unit vectors either as std::array<double,3> or as af::array.
    /// Input is normalized in ctor of UnitVectorOrArray class
    const UnitVectorOrArray c1, c2, c3;

    af::array h(const State& state);
    double get_cpu_time() { return 0; } // TODO//!< accumulated heff computation time in [s]
    double E(const State& state);
    double E(const State& state, const af::array& h);

  private:
    std::array<af::array, 3> h_1to3(const State& state);
};
} // namespace magnumafcpp
