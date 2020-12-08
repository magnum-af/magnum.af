#pragma once
#include <math.h>

namespace magnumafcpp {

namespace constants {
// The following values are obtained from CODATA/NIST as of 04.03.2019:

constexpr double mu0 = 4e-7 * M_PI; ///< [H/m] magnetic constant mu_0

constexpr double gamma = 1.760859644e11 * mu0; ///< [m A^-1 s^-1] gyromagnetic ratio gamma. Electron gyromagnetic
                                               ///< ratio gamma_e [s^-1 T^-1] times mu_0 [H/m], where [T/H] is [A/m^2]

constexpr double mu_b = 9.274009994e-24; ///< [J/T] Bohr magneton mu_bohr

constexpr double e_neg = -1.6021766208e-19; ///< [C] Elementary charge e (including minus sign)

constexpr double e_abs = 1.6021766208e-19; ///< [C] Absolute value of the elementary charge e

constexpr double kb = 1.38064852e-23; ///< [J/K] Boltzmann constant kb

constexpr double hbar = 1.0545718e-34; ///< [J s] Reduced Planck constant
} // namespace constants

namespace conversion {
/// Conversion of Ampere per meter [A/m] to Tesla [T]
template <typename T> T Am_to_Tesla(T Ampere_per_meter) { return Ampere_per_meter * constants::mu0; }

/// Conversion of Tesla [T] to Ampere per meter [A/m]
template <typename T> T Tesla_to_Am(T Tesla) { return Tesla / constants::mu0; }

/// Conversion of Joule [J] to electron volt [eV]
template <typename T> T J_to_eV(T Joule) { return Joule / constants::e_abs; }

/// Conversion of electron volt [eV] to Joule [J]
template <typename T> T eV_to_J(T eV) { return eV * constants::e_abs; }
} // namespace conversion

} // namespace magnumafcpp
