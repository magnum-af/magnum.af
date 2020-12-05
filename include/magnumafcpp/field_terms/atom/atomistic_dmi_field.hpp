#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"

namespace magnumafcpp {

class AtomisticDmiField : public AtomisticTermBase {
  public:
    AtomisticDmiField(const double D_atom, std::array<double, 3> D_atom_axis);
    AtomisticDmiField(const double D_atom, double D_atom_axis_x, double D_atom_axis_y, double D_atom_axis_z);
    virtual af::array h(const State& state) const override;

  private:
    const double D_atom;                     //!< Atomistic DMI energy [J]
    const std::array<double, 3> D_atom_axis; //!< Atomistic DMI axis
};
} // namespace magnumafcpp
