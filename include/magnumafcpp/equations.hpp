#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "util/func.hpp" // for cross4
namespace magnumafcpp::equations {
inline af::array LLG(double alpha, const af::array& m, const af::array& h_eff) {
    const af::array m_x_h = cross4(m, h_eff);
    return -constants::gamma / (1. + std::pow(alpha, 2)) * m_x_h -
           alpha * constants::gamma / (1. + std::pow(alpha, 2)) * cross4(m, m_x_h);
}

inline af::array dampingless_LLG(double alpha, const af::array& m, const af::array& h_eff) {
    return -alpha * constants::gamma / (1. + std::pow(alpha, 2)) * cross4(m, cross4(m, h_eff));
}
} // namespace magnumafcpp::equations
