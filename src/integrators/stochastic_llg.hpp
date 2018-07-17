#ifndef STOCHASTIC_LLG_H
#define STOCHASTIC_LLG_H
#include "arrayfire.h"
#include <memory>
#include "../llg_terms/LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
#include "stochastic_integrator.hpp"
using namespace af;

class Stochastic_LLG : public Stochastic_Integrator {
    public:
        Stochastic_LLG(State state, std::vector<std::shared_ptr<LLGTerm>> terms, const double d, std::string smode): Stochastic_Integrator (state, terms, d, smode){}
        double E(const State& state); //Energy calculation
    private:
        array fheff(const array& m);
        array detfdmdt(const array& m);//only for reference in detRK4
        array stochfdmdt(const array& m, const array& h_th);
};

#endif
