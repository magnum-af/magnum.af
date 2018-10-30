#include "new_llg.hpp"

NewLlg::NewLlg(std::string scheme, Controller controller, bool dissipation_term_only) : AdaptiveRungeKutta(scheme, controller), dissipation_term_only(dissipation_term_only) {
};

NewLlg::NewLlg(LlgTerms llgterms, std::string scheme, Controller controller, bool dissipation_term_only) : llgterms(llgterms), AdaptiveRungeKutta(scheme, controller), dissipation_term_only(dissipation_term_only) {
};

af::array NewLlg::fheff(const State& state){
  af::array solution = constant(0.,state.mesh.dims, f64);
  af::timer timer_heff = af::timer::start();

  for(unsigned i=0;i<llgterms.size();++i){
    solution+=llgterms[i]->h(state);
  }
  time_heff+=af::timer::stop(timer_heff);
  return solution;
}

af::array NewLlg::f(const State& state){
  //calls_fdmdt++;
  //timer_fdmdt=timer::start();
  if(dissipation_term_only){
    af::array heff=fheff(state);
    return - state.param.alpha*state.param.gamma/(1.+pow(state.param.alpha,2)) * cross4(state.m, cross4(state.m, heff));
  }
  else{
    af::array heff=fheff(state);
    return - state.param.gamma/(1.+pow(state.param.alpha,2)) * cross4(state.m, heff) - state.param.alpha*state.param.gamma/(1.+pow(state.param.alpha,2)) * cross4(state.m, cross4(state.m, heff));
  }
  //time_fdmdt+=af::timer::stop(timer_fdmdt);
}

//Energy calculation
double NewLlg::E(const State& state){
  double solution = 0.;
  for(unsigned i=0;i<llgterms.size();++i){
    solution+=llgterms[i]->E(state);
  }
  return solution;
}

void NewLlg::relax(State& state, const double precision, const int iloop, const int iwritecout){
    af::timer t = af::timer::start();
    double E_prev=1e20;
    while (fabs((E_prev-E(state))/E_prev) > precision){
        E_prev=E(state);
        for ( int i = 0; i<iloop; i++){
            step(state);
        }
        if( state.steps % iwritecout == 0) std::cout << "step " << state.steps << " rdiff= " << fabs((E_prev - E(state))/E_prev) << std::endl;
    }
    std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) << ", current llg steps = " << state.steps << std::endl; 
}

long int NewLlg::get_fheff_addr(const State& state){
    //std::vector<af::array> get_fheff_addr_temp_array;
    //get_fheff_addr_temp_array.push_back(fheff(state));
    //get_fheff_addr_temp_array.back().lock(); 
    //return (long int) get_fheff_addr_temp_array.back().get();
    //TODO elaborate best way to handle 
    fheff_tmp=fheff(state);
    return (long int) fheff_tmp.get();
}
