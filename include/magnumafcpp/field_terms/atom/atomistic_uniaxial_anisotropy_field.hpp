#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"
#include <array>

namespace magnumafcpp {

class AtomisticUniaxialAnisotropyField : public AtomisticTermBase {
  public:
    AtomisticUniaxialAnisotropyField(const double K_atom, std::array<double, 3> K_atom_axis = {0, 0, 1});
    AtomisticUniaxialAnisotropyField(const double K_atom, double K_atom_axis_x, double K_atom_axis_y,
                                     double K_atom_axis_z);
    af::array h(const State& state) const override;

  private:
    double K_atom;                     //!< Atomistic anisotropy energy in [J]
    std::array<double, 3> K_atom_axis; //!< Atomistic anisotropy axis
    std::array<double, 3> get_normalized_vector(std::array<double, 3> vector);
};
} // namespace magnumafcpp
