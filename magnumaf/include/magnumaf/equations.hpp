#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "math.hpp"
namespace magnumaf::equations {

// Warning: nasty bug: using pow2 in this scope as name for template will prefer call to af::pow2 over tempalte named
// pow2!!! This is due to the fact that overload resution comes before template specializations, thus concrete functions
// have precedence over templated ones. Therefore either introduce sub-namespace util::pow2, or choose other name, best
// would be both. However, function not have to be templated at all for this case as we need specialization, better
// pow_2 is overloaded function anyway, then 'call is ambiguous' error kicks in, preventing this bug.
// Takeaway: overload if possible, spezialization only if needed.
namespace detail {
inline double pow_2(const double& alpha) { return std::pow(alpha, 2.0); }
inline af::array pow_2(const af::array& alpha) { return af::pow(alpha, 2.0); }
} // namespace detail

/// LLG precession term
template <typename T> af::array LLG_precession(T alpha, const af::array& m_x_h) {
    return -constants::gamma / (1. + detail::pow_2(alpha)) * m_x_h;
}

/// LLG damping term
template <typename T> af::array LLG_damping(T alpha, const af::array& m, const af::array& m_x_h) {
    // TODO NOTE:: using pow2 as name for above template calls af::pow2!!!
    return -alpha * constants::gamma / (1. + detail::pow_2(alpha)) * math::cross4(m, m_x_h);
}

/// LLG equation
template <typename T> af::array LLG(T alpha, const af::array& m, const af::array& h_eff) {
    const af::array m_x_h = math::cross4(m, h_eff);
    return LLG_precession(alpha, m_x_h) + LLG_damping(alpha, m, m_x_h);
}

} // namespace magnumaf::equations
