#include "new_llg.hpp"

LLGIntegrator::LLGIntegrator(double alpha, std::string scheme, Controller controller, bool dissipation_term_only) : alpha(alpha), AdaptiveRungeKutta(scheme, controller), dissipation_term_only(dissipation_term_only) {
};

LLGIntegrator::LLGIntegrator(double alpha, LlgTerms llgterms, std::string scheme, Controller controller, bool dissipation_term_only) : alpha(alpha), AdaptiveRungeKutta(scheme, controller), llgterms(llgterms),  dissipation_term_only(dissipation_term_only) {
};

af::array LLGIntegrator::fheff(const State& state){
  af::array solution = constant(0.,state.mesh.dims, f64);
  af::timer timer_heff = af::timer::start();

  for(unsigned i=0;i<llgterms.size();++i){
    solution+=llgterms[i]->h(state);
  }
  time_heff+=af::timer::stop(timer_heff);
  return solution;
}

af::array LLGIntegrator::f(const State& state){
  //calls_fdmdt++;
  //timer_fdmdt=timer::start();
  if(dissipation_term_only){
    af::array heff=fheff(state);
    return - alpha*constants::gamma/(1.+pow(alpha,2)) * cross4(state.m, cross4(state.m, heff));
  }
  else{
    af::array heff=fheff(state);
    return - constants::gamma/(1.+pow(alpha,2)) * cross4(state.m, heff) - alpha*constants::gamma/(1.+pow(alpha,2)) * cross4(state.m, cross4(state.m, heff));
  }
  //time_fdmdt+=af::timer::stop(timer_fdmdt);
}

//Energy calculation
double LLGIntegrator::E(const State& state){
  double solution = 0.;
  for(unsigned i=0;i<llgterms.size();++i){
    solution+=llgterms[i]->E(state);
  }
  return solution;
}

void LLGIntegrator::relax(State& state, const double precision, const int iloop, const int iwritecout){
    af::timer t = af::timer::start();
    double E_prev=1e20;
    while (fabs((E_prev-E(state))/E_prev) > precision){
        E_prev=E(state);
        for ( int i = 0; i<iloop; i++){
            step(state);
        }
        if( iwritecout > 0 and state.steps % iwritecout == 0) printf("LLGIntegrator: Relax: step %llu, rdiff= %e", state.steps, fabs((E_prev - E(state))/E_prev));
    }
    printf("timerelax [af-s]: %e . Current state.steps= %llu and state.t = %e", state.t, state.steps, af::timer::stop(t)); 
}

long int LLGIntegrator::get_fheff_addr(const State& state){
    //std::vector<af::array> get_fheff_addr_temp_array;
    //get_fheff_addr_temp_array.push_back(fheff(state));
    //get_fheff_addr_temp_array.back().lock(); 
    //return (long int) get_fheff_addr_temp_array.back().get();
    
    //TODO tempfix for wrapping, elaborate other solution
    fheff_tmp=fheff(state);
    return (long int) fheff_tmp.get();
    //return (long int) fheff(state).get();
}
