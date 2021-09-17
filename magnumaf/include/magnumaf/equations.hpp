#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "math.hpp"
namespace magnumaf::equations {

template <typename T> T pow2(const T& alpha);
template <> inline double pow2<double>(const double& alpha){
    return std::pow(alpha, 2.0);
}
template <> inline af::array pow2<af::array>(const af::array& alpha){
    return af::pow(alpha, 2.0);
}

/// LLG precession term
template <typename T>
af::array LLG_precession(T alpha, const af::array& m_x_h) {
    return -constants::gamma / (1. + pow2(alpha)) * m_x_h;
}

/// LLG damping term
template <typename T>
af::array LLG_damping(T alpha, const af::array& m, const af::array& m_x_h) {
    return -alpha * constants::gamma / (1. + pow2(alpha)) * math::cross4(m, m_x_h);
}

/// LLG equation
template <typename T>
af::array LLG(T alpha, const af::array& m, const af::array& h_eff) {
    const af::array m_x_h = math::cross4(m, h_eff);
    return LLG_precession(alpha, m_x_h) + LLG_damping(alpha, m, m_x_h);
}

} // namespace magnumaf::equations
