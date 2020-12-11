#pragma once
#include "arrayfire.h"
#include <functional>
#include <utility>

namespace magnumafcpp {

/// RK4 Method
/// @param t0 Initial time
/// @param y0 Initial config
/// @param dt stepsize
/// @param f  Callback function f(t,y,...)=dy/dt
/// @param args Variadic args ... for f(t, y, ... )
template <typename D, typename T, typename F, typename... Args>
// std::pair<D, T> RK4(D t0, T y0, D dt, F f, Args... args) {
std::pair<D, T> RK4(const D t0, const D dt, const T& y0, const F& f, const Args&... args) {
    // Running k, ksum
    // k1
    T k = f(t0, y0, args...);
    T ksum = k;
    // k2
    k = f(t0 + dt / 2., y0 + dt * k / 2., args...);
    ksum += 2 * k;
    // k3
    k = f(t0 + dt / 2., y0 + dt * k / 2., args...);
    ksum += 2 * k;
    // k4
    ksum += f(t0 + dt, y0 + dt * k, args...);
    return {t0 + dt, y0 + 1. / 6. * dt * (ksum)};
}
} // namespace magnumafcpp
