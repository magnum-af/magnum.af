#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"
#include <array>

namespace magnumafcpp
{

class AtomisticUniaxialAnisotropyField : public LLGTerm
{
public:
    //Field contribution
    af::array h(const State &state);
    //CPU time
    double get_cpu_time() { return cpu_time; }

    AtomisticUniaxialAnisotropyField(const double K_atom, std::array<double, 3> K_atom_axis = {0, 0, 1});

    af::array eu; //Uniaxial anisotropy normal vector

    double cpu_time{0.};
    const double K_atom;
    const std::array<double, 3> K_atom_axis; //!< Anisotropy axis
private:
    std::array<double, 3> get_normalized_vector(std::array<double, 3> vector);
};
} // namespace magnumafcpp
