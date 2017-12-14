#include "stochastic_llg.hpp"
Stochastic_LLG::Stochastic_LLG (State in, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in) : Fieldterms(Fieldterms_in),  param(in.param), mesh(in.mesh){
}

array Stochastic_LLG::fheff(const array& m){
    State temp(mesh,param,m);
    array solution = Fieldterms[0]->h(temp);
    for(unsigned i=1;i<Fieldterms.size();++i){
        solution+=Fieldterms[i]->h(temp);
    }
    return solution;
}

//array Stochastic_LLG::fdmdt(const array& m){
//    array heff = fheff(m);
//    array cross_temp = cross4(m, heff);
//    return - param.gamma/(1.+pow(param.alpha,2)) * cross_temp - param.alpha*param.gamma/(1.+pow(param.alpha,2)) * cross4(m, cross_temp);
//}


//array Stochastic_LLG::rk4(const array& m, const double dt)
//{
//  array k1   =  dt * fdmdt(m                               );
//  array k2   =  dt * fdmdt(m + 1./2.*k1                    );
//  array k3   =  dt * fdmdt(m            + 1./2.*k2         );
//  array k4   =  dt * fdmdt(m                       +    k3 );
//
//  return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
//}

array Stochastic_LLG::fdmdt(const array& m, const array& heff){
    array cross_temp = cross4(m, heff);
    return - param.gamma/(1.+pow(param.alpha,2)) * cross4(m, heff) - param.alpha*param.gamma/(1.+pow(param.alpha,2)) * cross4(m, cross4(m, heff));
}

array Stochastic_LLG::rk4(const array& m, const double dt)
{
  array heff =       fheff(m                                              );
  array k1   =  dt * fdmdt(m                                        , heff);
  heff =       fheff(m + 1./2.*k1                                   );
  array k2   =  dt * fdmdt(m + 1./2.*k1                             , heff);
  heff =       fheff(m            + 1./2.*k2                        );
  array k3   =  dt * fdmdt(m            + 1./2.*k2                  , heff);
  heff =       fheff(m                       +    k3                );
  array k4   =  dt * fdmdt(m                       +    k3          , heff);

  return            (          k1 +    2.*k2 + 2.*k3 + k4    ) / 6.;
}

//array Stochastic_LLG::StemiImplicitHeun(const array& m, const double dt){
//    if(calls > 0) { 
//        array m_it = m + ( m - m_prev)/2.;
//    }
//    else {
//        m_prev = m; 
//        std::cout<<"CALLS !>0"<<std::endl;
//    }
//    array heff = fheff(m);
//    for (int i = 0; i < 5; i++){
//       m_it = m + fdmdt(m_it,heff) * dt/2.;
//    }
//    m_prev = m + fdmdt(m_it,heff) * dt;
//    return m_prev;
//}

void Stochastic_LLG::step(State& state, const double dt){
    state.m += rk4(state.m,dt);
    //state.m = StemiImplicitHeun(state.m,dt);
    state.m = renormalize(state.m);
    state.t+=dt;
    calls ++;
}
