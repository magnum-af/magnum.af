#pragma once
#include <tuple>
namespace magnumafcpp {

///
/// /param
template <typename D, typename T, typename F, typename... Args>
std::tuple<D, T, T, T> BS23(const D t_n, const D dt, const T& y_n, const T& k1, const F& f, const Args&... args) {
    const T k2 = dt * f(t_n + 1. / 2. * dt, y_n + 1. / 2. * k1, args...);
    const T k3 = dt * f(t_n + 3. / 4. * dt, y_n + 3. / 4. * k2, args...);

    const T dy = 2. / 9. * k1 + 1. / 3. * k2 + 4. / 9. * k3;
    const T k4 = dt * f(t_n + dt, y_n + dy, args...);
    const T rk_error = dy - (7. / 24. * k1 + 1. / 4. * k2 + 1. / 3. * k3 + 1. / 8. * k4);
    return {t_n + dt, y_n + dy, rk_error, k4};
}

template <typename D, typename T, typename F, typename... Args>
std::tuple<D, T, T, T> BS23(const D t_n, const D dt, const T& y_n, const F& f, const Args&... args) {
    const T k1 = dt * f(t_n, y_n, args...);
    return {BS23(t_n, dt, y_n, k1, f, args...)};
}

} // namespace magnumafcpp
