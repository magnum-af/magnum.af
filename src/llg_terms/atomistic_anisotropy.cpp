#include "atomistic_anisotropy.hpp"
#include "../func.hpp"

namespace magnumafcpp {

std::array<double, 3> AtomisticUniaxialAnisotropyField::get_normalized_vector(std::array<double, 3> vector) {
    double norm = sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
    return std::array<double, 3>{vector[0] / norm, vector[1] / norm, vector[2] / norm};
}

// Ref Master Thesis Stifano
// eq (19)
// Han=N*K/2 - K/2 Sum_i(m_i*ez)^2

AtomisticUniaxialAnisotropyField::AtomisticUniaxialAnisotropyField(const double K_atom,
                                                                   std::array<double, 3> K_atom_axis)
    : K_atom(K_atom), K_atom_axis(get_normalized_vector(K_atom_axis)) {}

AtomisticUniaxialAnisotropyField::AtomisticUniaxialAnisotropyField(const double K_atom, double K_atom_axis_x,
                                                                   double K_atom_axis_y, double K_atom_axis_z)
    : K_atom(K_atom),
      K_atom_axis(get_normalized_vector(std::array<double, 3>{K_atom_axis_x, K_atom_axis_y, K_atom_axis_z})) {}

af::array AtomisticUniaxialAnisotropyField::h(const State& state) {
    af::timer timer_anisotropy = af::timer::start();
    // Normal vector
    eu = af::array(state.m.dims(), f64);
    eu(af::span, af::span, af::span, 0) = K_atom_axis[0];
    eu(af::span, af::span, af::span, 1) = K_atom_axis[1];
    eu(af::span, af::span, af::span, 2) = K_atom_axis[2];
    af::array anisotropy = eu * state.m;
    anisotropy = af::sum(anisotropy, 3);
    anisotropy = af::tile(anisotropy, 1, 1, 1, 3);

    if (state.afsync) {
        af::sync();
    }
    cpu_time += af::timer::stop(timer_anisotropy);
    return 2 * K_atom / (constants::mu0 * state.Ms) * anisotropy * eu;
}

} // namespace magnumafcpp
