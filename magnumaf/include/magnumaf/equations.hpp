#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "util/double_or_array.hpp"
#include "math.hpp"
namespace magnumaf::equations {

/// LLG precession term
inline af::array LLG_precession(util::DoubleOrArray const& alpha, const af::array& m_x_h) {
    const auto alpha_pow2 = pow2_vec(alpha, m_x_h.dims(), m_x_h.type());
    return -constants::gamma / (1. + alpha_pow2) * m_x_h;
}

/// LLG damping term
inline af::array LLG_damping(util::DoubleOrArray const& alpha, const af::array& m, const af::array& m_x_h) {
    const auto alpha_pow2 = pow2_vec(alpha, m_x_h.dims(), m_x_h.type());
    return alpha * (-constants::gamma / (1. + alpha_pow2 * math::cross4(m, m_x_h)));
}

/// LLG equation
inline af::array LLG(util::DoubleOrArray const& alpha, const af::array& m, const af::array& h_eff) {
    // TODO //  alpha.print();
    const af::array m_x_h = math::cross4(m, h_eff);
    return LLG_precession(alpha, m_x_h) + LLG_damping(alpha, m, m_x_h);
}

} // namespace magnumaf::equations
