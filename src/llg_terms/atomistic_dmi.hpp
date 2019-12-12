#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumafcpp{

class AtomisticDmiField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    //AtomisticDmiField ();
    AtomisticDmiField (const double D_atom, std::array<double, 3> D_atom_axis);

    double     cpu_time{0.};
    const double D_atom;  //!< Atomistic DMI energy [J]
    const std::array<double, 3> D_atom_axis; //!< Atomistic DMI axis
};
}// namespace magnumafcpp
