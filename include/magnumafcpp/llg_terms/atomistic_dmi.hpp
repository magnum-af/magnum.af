#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "atomistic_term_base.hpp"

namespace magnumafcpp {

class AtomisticDmiField : public AtomisticTermBase {
  public:
    // Field contribution
    af::array h(const State& state);
    // CPU time
    double get_cpu_time() { return cpu_time; }

    // AtomisticDmiField ();
    AtomisticDmiField(const double D_atom, std::array<double, 3> D_atom_axis);
    AtomisticDmiField(const double D_atom, double D_atom_axis_x, double D_atom_axis_y, double D_atom_axis_z);

    double cpu_time{0.};
    const double D_atom;                     //!< Atomistic DMI energy [J]
    const std::array<double, 3> D_atom_axis; //!< Atomistic DMI axis
};
} // namespace magnumafcpp
