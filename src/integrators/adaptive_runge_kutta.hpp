#ifndef ADAPTIVE_INTEGRATOR_H
#define ADAPTIVE_INTEGRATOR_H
#include "arrayfire.h"
#include "controller.hpp"
#include "../func.hpp"
#include "../state.hpp"

typedef State (*callback_function)(const State& state);

class AdaptiveRungeKutta {
    public:
        AdaptiveRungeKutta(callback_function f, std::string scheme = "RKF45", Controller controller = Controller());
        void step(af::array& m, double& t);
    private:
        callback_function f;
        const std::string scheme; //Integration scheme s.a. RKF45, DP45, ...
        Controller controller;
        af::array RKF45(const af::array& m, const double t, const double dt, double& err);
        af::array DP45(const af::array& m, const double t, const double dt, double& err);
        af::array BS45(const af::array& m, const double t, const double dt , double& err);
        double h{1.01e-15}; //step size of RK used in controller
        double err{0}; // error for stepsize controller
};

#endif
