#include "stochastic_llg.hpp"
#include "../llg_terms/LLGTerm.hpp"
#include "../func.hpp"
#include <memory>

namespace magnumafcpp{


//Energy calculation
double Stochastic_LLG::E(const State& state){
    double solution = 0.;
    for(unsigned i=0;i<Fieldterms.size();++i){
        solution+=Fieldterms[i]->E(state);
    }
    return solution;
}

af::array Stochastic_LLG::fheff(const State& state){
    af::array solution = Fieldterms[0]->h(state);
    for(unsigned i=1;i<Fieldterms.size();++i){
        solution+=Fieldterms[i]->h(state);
    }
    return solution;
}

af::array Stochastic_LLG::detfdmdt(const State& state){
    fdmdt_calls++;
    const af::array heff = fheff(state);
    const af::array cross_temp = cross4(state.m, heff);
    return - constants::gamma/(1.+pow(this->alpha, 2)) * cross_temp - this->alpha*constants::gamma/(1.+pow(this->alpha, 2)) * cross4(state.m, cross_temp);
}

af::array Stochastic_LLG::stochfdmdt(const State& state, const af::array& h_th){
    stochfdmdt_calls++;
    const af::array h = fheff(state) + h_th;
    const af::array cross_temp = cross4(state.m, h);
    return  - constants::gamma/(1.+pow(this->alpha, 2)) * cross_temp - this->alpha*constants::gamma/(1.+pow(this->alpha, 2)) * cross4(state.m, cross_temp);
}

}// namespace magnumafcpp
