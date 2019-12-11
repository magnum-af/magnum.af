#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumafcpp
{

class AtomisticUniaxialAnisotropyField : public LLGTerm
{
public:
    //Field contribution
    af::array h(const State &state);
    //Energy contribution
    double E(const State &state);
    double E(const State &state, const af::array &h); ///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time() { return cpu_time; }

    AtomisticUniaxialAnisotropyField(const double K_atom, std::array<double, 3> K_atom_axis = {0, 0, 1});

    af::array eu; //Uniaxial anisotropy normal vector

    double cpu_time{0.};
    const double K_atom;
    const std::array<double, 3> K_atom_axis; //!< Anisotropy axis
};
} // namespace magnumafcpp
