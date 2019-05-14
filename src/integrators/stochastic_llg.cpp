#include "stochastic_llg.hpp"
using namespace af;

//Energy calculation
double Stochastic_LLG::E(const State& state){
    double solution = 0.;
    for(unsigned i=0;i<Fieldterms.size();++i){
        solution+=Fieldterms[i]->E(state);
    }
    return solution;
}

array Stochastic_LLG::fheff(const array& m){
    State temp(mesh,material,m);
    array solution = Fieldterms[0]->h(temp);
    for(unsigned i=1;i<Fieldterms.size();++i){
        solution+=Fieldterms[i]->h(temp);
    }
    return solution;
}

array Stochastic_LLG::detfdmdt(const array& m){
    fdmdt_calls++;
    const array heff = fheff(m);
    const array cross_temp = cross4(m, heff);
    return - constants::gamma/(1.+pow(material.alpha,2)) * cross_temp - material.alpha*constants::gamma/(1.+pow(material.alpha,2)) * cross4(m, cross_temp);
}

array Stochastic_LLG::stochfdmdt(const array& m, const array& h_th){
    stochfdmdt_calls++;
    const array h = fheff(m) + h_th;
    const array cross_temp = cross4(m, h);
    return  - constants::gamma/(1.+pow(material.alpha,2)) * cross_temp - material.alpha*constants::gamma/(1.+pow(material.alpha,2)) * cross4(m, cross_temp);
}
 
