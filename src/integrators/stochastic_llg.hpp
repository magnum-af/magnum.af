#ifndef STOCHASTIC_LLG_H
#define STOCHASTIC_LLG_H
#include "arrayfire.h"
#include <memory>
#include "../llg_terms/LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
#include "stochastic_integrator.hpp"

class Stochastic_LLG : public Stochastic_Integrator {
    public:
        Stochastic_LLG(double alpha, double T, double dt, State state, std::vector<std::shared_ptr<LLGTerm>> terms, std::string smode): Stochastic_Integrator (alpha, T, dt, state, terms, smode){}
        double E(const State& state); //Energy calculation
    private:
        af::array fheff(const State&);
        af::array detfdmdt(const State&);//only for reference in detRK4
        af::array stochfdmdt(const State&, const af::array& h_th);
};

#endif
