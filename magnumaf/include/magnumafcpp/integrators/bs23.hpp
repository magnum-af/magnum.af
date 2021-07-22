#pragma once
#include <tuple>

/// Runge-Kutta methods for the integration of ordinary differential equations.
namespace magnumaf::integrator {

/// Bogacki–Shampine 2nd/3rd order method using FSAL property
/// \param k1 FSAL k4 stage of previous step, used as k1
template <typename D, typename T, typename F, typename... Args>
std::tuple<D, T, T, T> BS23(const D t_n, const D dt, const T& y_n, const T& k1, const F& f, const Args&... args) {
    const T k2 = dt * f(t_n + 1. / 2. * dt, y_n + 1. / 2. * k1, args...);
    const T k3 = dt * f(t_n + 3. / 4. * dt, y_n + 3. / 4. * k2, args...);

    const T dy = 2. / 9. * k1 + 1. / 3. * k2 + 4. / 9. * k3; // 3rd order estimate (diff, w.o. y_n)
    const T k4 = dt * f(t_n + dt, y_n + dy, args...);
    const T dz = 7. / 24. * k1 + 1. / 4. * k2 + 1. / 3. * k3 + 1. / 8. * k4; // 2nd order estimate (diff, w.o. y_n)
    const T rk_error = dy - dz;
    return {t_n + dt, y_n + dy, rk_error, k4};
}

/// Bogacki–Shampine 2nd/3rd order method.
/// Uses 3rd order (dy) for propagation and 2nd order (rk_error) for error estimate.
/// returns {t_n + dt, y_n + dy, rk_error, k4}
template <typename D, typename T, typename F, typename... Args>
std::tuple<D, T, T, T> BS23(const D t_n, const D dt, const T& y_n, const F& f, const Args&... args) {
    const T k1 = dt * f(t_n, y_n, args...);
    return BS23(t_n, dt, y_n, k1, f, args...);
}

} // namespace magnumaf::integrator
