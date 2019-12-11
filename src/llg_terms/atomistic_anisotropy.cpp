#include "atomistic_anisotropy.hpp"
#include "../func.hpp"
#include "../util/util.hpp"

namespace magnumafcpp
{

//Ref Master Thesis Stifano
//eq (19)
//Han=N*K/2 - K/2 Sum_i(m_i*ez)^2

AtomisticUniaxialAnisotropyField::AtomisticUniaxialAnisotropyField(const double K_atom, std::array<double, 3> K_atom_axis) : K_atom(K_atom), K_atom_axis(get_normalized_vector(K_atom_axis))
{
}

af::array AtomisticUniaxialAnisotropyField::h(const State &state)
{
    af::timer timer_anisotropy = af::timer::start();
    //Normal vector
    eu = af::array(state.mesh.n0, state.mesh.n1, state.mesh.n2, 3, f64);
    eu(af::span, af::span, af::span, 0) = K_atom_axis[0];
    eu(af::span, af::span, af::span, 1) = K_atom_axis[1];
    eu(af::span, af::span, af::span, 2) = K_atom_axis[2];
    af::array anisotropy = eu * state.m;
    anisotropy = af::sum(anisotropy, 3);
    anisotropy = af::tile(anisotropy, 1, 1, 1, 3);

    if (state.afsync)
    {
        af::sync();
    }
    cpu_time += af::timer::stop(timer_anisotropy);
    return 2 * K_atom / (constants::mu0 * state.Ms) * anisotropy * eu;
}

}// namespace mangumafcpp
