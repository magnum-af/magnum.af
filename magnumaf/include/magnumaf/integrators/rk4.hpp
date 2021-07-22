#pragma once
#include "arrayfire.h"
#include <utility>

namespace magnumaf {

/// RK4 Method
/// Returns {tn + dt, yn + dy}
/// @param tn Initial time
/// @param yn Initial config
/// @param dt stepsize
/// @param f  Callback function f(t,y,...)=dy/dt
/// @param args Variadic args ... for f(t, y, ... )
// Assuming D is mostly double and hence cheap to copy we take it by value
template <typename D, typename T, typename F, typename... Args>
std::pair<D, T> RK4(const D tn, const D dt, const T& yn, const F& f, const Args&... args) {
    // k1
    T k = f(tn, yn, args...);
    T ksum = k;
    // k2
    k = f(tn + dt / 2., yn + dt * k / 2., args...);
    ksum += 2 * k;
    // k3
    k = f(tn + dt / 2., yn + dt * k / 2., args...);
    ksum += 2 * k;
    // k4
    ksum += f(tn + dt, yn + dt * k, args...);
    return {tn + dt, yn + 1. / 6. * dt * (ksum)};
}
} // namespace magnumaf
