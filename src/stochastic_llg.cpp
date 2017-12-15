#include "stochastic_llg.hpp"
Stochastic_LLG::Stochastic_LLG (State in, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in) : Fieldterms(Fieldterms_in),  param(in.param), mesh(in.mesh), m_prev(in.m){
}

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
    array heff = fheff(m);
    array cross_temp = cross4(m, heff);
    return - param.gamma/(1.+pow(param.alpha,2)) * cross_temp - param.alpha*param.gamma/(1.+pow(param.alpha,2)) * cross4(m, cross_temp);
}

double Stochastic_LLG::cpu_time(){
  double cpu_time = 0.;
  for(unsigned i=0;i<Fieldterms.size();++i){
    cpu_time+=Fieldterms[i]->get_cpu_time();
  }
  return cpu_time;
}


array Stochastic_LLG::rk4(const array& m, const double dt)
{
  array k1   =  dt * fdmdt(m                               );
  array k2   =  dt * fdmdt(m + 1./2.*k1                    );
  array k3   =  dt * fdmdt(m            + 1./2.*k2         );
  array k4   =  dt * fdmdt(m                       +    k3 );

  return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
}

//array Stochastic_LLG::StemiImplicitHeun(const array& m, const double dt){
//    timer_stoch = timer::start();
//    array m_it = m + ( m - m_prev)/2.;
//    for (int i = 0; i < 5; i++){
//       m_it = m + fdmdt(m_it) * dt/2.;
//    }
//    m_prev = m;
//    std::cout<<" CPU TIME  = "<<cpu_time()<<"\n"<<std::endl;
//    return m + fdmdt(m_it) * dt;
//}

array Stochastic_LLG::StemiImplicitHeun(const array& m, const double dt){
    timer_stoch = timer::start();
    //Now in constructor
    //if(calls == 0) { 
    //    m_prev = m; 
    //    std::cout<<"CALLS == 0"<<std::endl;
    //}
    //TODO this is somehow super slow
    //
//    array m_it = m + ( m - m_prev)/2.;
//    for (int i = 0; i < 5; i++){
//       m_it = m + fdmdt(m_it) * dt/2.;
//    }
//    m_prev = m;
//    std::cout<<" CPU TIME  = "<<cpu_time()<<"\n"<<std::endl;
//    //return m + fdmdt(m_it) * dt;
//    return  fdmdt(m_it) * dt;
    
    array m0= 1.5* m - 0.5 * m_prev;
    m_prev = m;
    array m1= m + dt/2. * fdmdt (m0);
    array m2= m + dt/2. * fdmdt (m1);
    array m3= m + dt/2. * fdmdt (m2);
    array m4= m + dt/2. * fdmdt (m3);
    array m5= m + dt/2. * fdmdt (m4);
    array result = dt * fdmdt(m5);
    time += timer::stop(timer_stoch);
    std::cout<<" CPU TIME  = "<<cpu_time()<<"\n"<<std::endl;
    //return dt *fdmdt(m5);
    return dt *fdmdt(m5);//TODO extremely slow for first ~100 runs, then fast

//    array k1   =  dt * fdmdt(m                               );
//    array k2   =  dt * fdmdt(m + 1./2.*k1                    );
//    array k3   =  dt * fdmdt(m            + 1./2.*k2         );
//    array k4   =  dt * fdmdt(m                       +    k3 );
//  
//    time += timer::stop(timer_stoch);
//    return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
}

void Stochastic_LLG::step(State& state, const double dt){
    //state.m += rk4(state.m,dt);
    state.m += StemiImplicitHeun(state.m,dt);
    state.m = renormalize(state.m);
    state.t+=dt;
    calls ++;
}


//array Stochastic_LLG::fdmdt(const array& m, const array& heff){
//    array cross_temp = cross4(m, heff);
//    return - param.gamma/(1.+pow(param.alpha,2)) * cross4(m, heff) - param.alpha*param.gamma/(1.+pow(param.alpha,2)) * cross4(m, cross4(m, heff));
//}
//
//array Stochastic_LLG::rk4(const array& m, const double dt)
//{
//  array heff =       fheff(m                                              );
//  array k1   =  dt * fdmdt(m                                        , heff);
//  heff =       fheff(m + 1./2.*k1                                   );
//  array k2   =  dt * fdmdt(m + 1./2.*k1                             , heff);
//  heff =       fheff(m            + 1./2.*k2                        );
//  array k3   =  dt * fdmdt(m            + 1./2.*k2                  , heff);
//  heff =       fheff(m                       +    k3                );
//  array k4   =  dt * fdmdt(m                       +    k3          , heff);
//
//  return            (          k1 +    2.*k2 + 2.*k3 + k4    ) / 6.;
//}

