#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "double_or_array.hpp"
#include "field_terms/micro/micro_term.hpp"
#include "unit_vector_or_array.hpp"
#include <array>

namespace magnumafcpp {

class CubicAnisotropyField : public MicroTerm {
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


    virtual af::array h(const State& state) const override;
    using MicroTerm::E;
    virtual double E(const State& state, const af::array& h) const override;

    DoubleOrArray Kc1, Kc2, Kc3; // First, second and third order cubic anisotropy constants in [J/m^3]
    /// Pairwise orthogonal unit vectors either as std::array<double,3> or as af::array.
    /// Input is normalized in ctor of UnitVectorOrArray class
    UnitVectorOrArray c1, c2, c3;

  private:
    std::array<af::array, 3> h_1to3(const State& state) const;
};
} // namespace magnumafcpp
