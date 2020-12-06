#pragma once
#include "arrayfire.h"
#include "field_terms/atom/atom_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

class AtomisticDmiField : public AtomTerm {
  public:
    AtomisticDmiField(double D_atom, std::array<double, 3> D_atom_axis);
    AtomisticDmiField(double D_atom, double D_atom_axis_x, double D_atom_axis_y, double D_atom_axis_z);
    virtual af::array h(const State& state) const override;

  private:
    double D_atom;                     //!< Atomistic DMI energy [J]
    std::array<double, 3> D_atom_axis; //!< Atomistic DMI axis
};
} // namespace magnumafcpp
