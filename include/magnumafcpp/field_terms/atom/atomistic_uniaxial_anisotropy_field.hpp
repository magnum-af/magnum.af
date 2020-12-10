#pragma once
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"
#include "state.hpp"
#include <array>

namespace magnumafcpp {

class AtomisticUniaxialAnisotropyField : public AtomTerm {
  public:
    AtomisticUniaxialAnisotropyField(const double K_atom, std::array<double, 3> K_atom_axis = {0, 0, 1});
    AtomisticUniaxialAnisotropyField(const double K_atom, double K_atom_axis_x, double K_atom_axis_y,
                                     double K_atom_axis_z);

  private:
    double K_atom;                     //!< Atomistic anisotropy energy in [J]
    std::array<double, 3> K_atom_axis; //!< Atomistic anisotropy axis
    std::array<double, 3> get_normalized_vector(std::array<double, 3> vector);

    virtual af::array impl_H_in_Apm(const State& state) const override;
};
} // namespace magnumafcpp
