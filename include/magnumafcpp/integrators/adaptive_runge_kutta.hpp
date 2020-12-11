#pragma once
#include "arrayfire.h"
#include "controller.hpp"
#include "state.hpp"
#include <stdio.h>
#include <string>

namespace magnumafcpp {

class AdaptiveRungeKutta {
  public:
    AdaptiveRungeKutta(std::string scheme_ = "RKF45", Controller controller_ = Controller(),
                       const bool normalize_ = true, const bool verbose = false);
    void step(State&);
    double get_time_allsteps() const { return time_allsteps_; }
    unsigned long long accumulated_steps{0}; //! accumulated integration steps, is incremented for each step of
                                             //! the integration scheme
    virtual ~AdaptiveRungeKutta(){};

  protected:
    af::array RKF45(const State& state, const double dt, double& err) const;

  private:
    virtual af::array f(const State& state) const = 0; // callback function s.a. LLG

    // settings params, would be const if not disabling move,copy ctors
    std::string scheme_; // could be const // Integration scheme s.a. RKF45, DP45, ...

    Controller controller_;

    bool normalize_; // could be const

    // Integration methods; not const specified do change k_FSAL
    af::array DP45(const State& state, const double dt, double& err);
    af::array BS45(const State& state, const double dt, double& err);
    af::array DP78(const State& state, const double dt, double& err) const;
    af::array BS23(const State& state, const double dt, double& err);

    double h_{1.01e-15}; // step size of RK used in controller
    double time_allsteps_{0};
    unsigned long long int step_calls_{0};

    af::array k_FSAL; // array which stores the last stage in methods with first
                      // same as last (FSAL) property
};

} // namespace magnumafcpp
