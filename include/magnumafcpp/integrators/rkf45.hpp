#pragma once
#include <tuple>

namespace magnumafcpp {

/// Runge-Kutta-Fehlberg 4th/5th order method
/// Returns {tn + dt, yn + dy(Ord5), rk_error ^= dy(Ord4)}
/// Using 5th order to propagate and 4th order for error estimate
// Assuming D is mostly double and hence cheap to copy we take it by value
template <typename D, typename T, typename F, typename... Args>
std::tuple<D, T, T> RKF45(const D tn, const D dt, const T& yn, const F& f, const Args&... args) {
    T k1 = dt * f(tn, yn, args...);
    T k2 = dt * f(tn + 1. / 4. * dt, yn + 1. / 4. * k1, args...);
    T k3 = dt * f(tn + 3. / 8. * dt, yn + 3. / 32. * k1 + 9 / 32. * k2, args...);
    T k4 = dt * f(tn + 12. / 13. * dt, yn + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3, args...);
    T k5 = dt * f(tn + dt, yn + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4, args...);
    T k6 = dt * f(tn + 1. / 2. * dt,
                  yn - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5, args...);
    T dy = 16. / 135. * k1 + 6656. / 12825. * k3 + 28561. / 56430. * k4 - 9. / 50. * k5 + 2. / 55. * k6;
    T rk_error = dy - (25. / 216. * k1 + 1408. / 2565. * k3 + 2197. / 4104. * k4 - 1. / 5. * k5);
    return {tn + dt, yn + dy, rk_error};
}
} // namespace magnumafcpp
