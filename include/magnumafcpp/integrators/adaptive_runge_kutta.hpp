#pragma once
#include "../state.hpp"
#include "controller.hpp"
#include "arrayfire.h"
#include <string>
#include <stdio.h>

namespace magnumafcpp
{

class AdaptiveRungeKutta
{
public:
    AdaptiveRungeKutta(std::string scheme_ = "RKF45", Controller controller_ = Controller(), const bool renormalize_ = true, const bool verbose = false);
    void step(State &);
    double get_time_allsteps() { return time_allsteps_; }
    virtual ~AdaptiveRungeKutta(){};

private:
    virtual af::array f(const State &state) = 0; // callback function s.a. LLG
    af::array RKF45(const State &state, const double dt, double &err);
    af::array DP45(const State &state, const double dt, double &err);
    af::array BS45(const State &state, const double dt, double &err);
    af::array DP78(const State &state, const double dt, double &err);
    af::array BS23(const State &state, const double dt, double &err);

    const std::string scheme_; //Integration scheme s.a. RKF45, DP45, ...
    Controller controller_;
    double h_{1.01e-15}; //step size of RK used in controller
    double err_{0};      // error for stepsize controller
    double time_allsteps_{0};
    unsigned long long int step_calls_{0};
    const bool renormalize_;

    af::array k_FSAL; // array which stores the last stage in methods with first same as last (FSAL) property
};

} // namespace magnumafcpp
