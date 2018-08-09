#ifndef ADAPTIVE_INTEGRATOR_H
#define ADAPTIVE_INTEGRATOR_H
#include "arrayfire.h"
#include "controller.hpp"
#include "../func.hpp"
#include "../state.hpp"

class AdaptiveRungeKutta {
    public:
        AdaptiveRungeKutta(std::string scheme = "RKF45", Controller controller = Controller());
        void step(State&);
    private:
        virtual af::array f(const State& state)=0; // callback function s.a. LLG
        const std::string scheme; //Integration scheme s.a. RKF45, DP45, ...
        Controller controller;
        af::array RKF45(const State& state, const double dt, double& err);
        af::array DP45(const af::array& m, const double t, const double dt, double& err);
        af::array BS45(const af::array& m, const double t, const double dt , double& err);
        double h{1.01e-15}; //step size of RK used in controller
        double err{0}; // error for stepsize controller
};

#endif
