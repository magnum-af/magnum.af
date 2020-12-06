#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "stochastic_integrator.hpp"

namespace magnumafcpp {

class Stochastic_LLG : public Stochastic_Integrator {
  public:
    Stochastic_LLG(double alpha, double T, double dt, State state, std::vector<std::unique_ptr<FieldTerm>> terms,
                   std::string smode)
        : Stochastic_Integrator(alpha, T, dt, state, std::move(terms), smode) {}
    virtual ~Stochastic_LLG() = default;
    double E(const State& state); // Energy calculation
  private:
    af::array fheff(const State&);
    af::array detfdmdt(const State&); // only for reference in detRK4
    af::array stochfdmdt(const State&, const af::array& h_th);
};

} // namespace magnumafcpp
