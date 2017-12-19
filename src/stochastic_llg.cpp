#include "stochastic_llg.hpp"
array Stochastic_LLG::fheff(const array& m){
    State temp(mesh,param,m);
    array solution = Fieldterms[0]->h(temp);
    for(unsigned i=1;i<Fieldterms.size();++i){
        solution+=Fieldterms[i]->h(temp);
    }
    return solution;
}

array Stochastic_LLG::fdmdt(const array& m){
    fdmdt_calls++;
    const array heff = fheff(m);
    const array cross_temp = cross4(m, heff);
    return - param.gamma/(1.+pow(param.alpha,2)) * cross_temp - param.alpha*param.gamma/(1.+pow(param.alpha,2)) * cross4(m, cross_temp);
}

array Stochastic_LLG::stochfdmdt(const array& m, const array& h_th){
    stochfdmdt_calls++;
    const array h = fheff(m) + h_th;
    const array cross_temp = cross4(m, h);
    return  - param.gamma/(1.+pow(param.alpha,2)) * cross_temp - param.alpha*param.gamma/(1.+pow(param.alpha,2)) * cross4(m, cross_temp);
}
 
