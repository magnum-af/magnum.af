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
        Stochastic_LLG(State state, std::vector<std::shared_ptr<LLGTerm>> terms, const double d, std::string smode): Stochastic_Integrator (state, terms, d, smode){}
        double E(const State& state); //Energy calculation
    private:
        af::array fheff(const af::array& m);
        af::array detfdmdt(const af::array& m);//only for reference in detRK4
        af::array stochfdmdt(const af::array& m, const af::array& h_th);
};

#endif
