#ifndef ADAPTIVE_INTEGRATOR_H
#define ADAPTIVE_INTEGRATOR_H
#include "arrayfire.h"
#include "controller.hpp"
#include "../func.hpp"

typedef af::array (*callback_function)(const af::array&);
class AdaptiveRungeKutta {
    public:
        AdaptiveRungeKutta(callback_function f, std::string scheme = "RKF45");
        af::array step(const af::array& m, const double dt);
    private:
        Controller controller = Controller();
        callback_function f;
        const std::string scheme; //Integration scheme s.a. RKF45, DP45, ...
        af::array RKF45(const af::array& m, const double dt, double& err);
};

#endif
