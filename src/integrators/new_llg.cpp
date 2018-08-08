#include "new_llg.hpp"

typedef array (*callback_function)(const State& state);
NewLlg::NewLlg(std::string scheme){
	State state(Mesh(0,0,0,0,0,0),Param(),af::constant(0,1));//TODO
    integrator =  AdaptiveRungeKutta(GetFPointer(), scheme);
};

af::array NewLlg::fheff(const State& state){
  af::array solution = constant(0.,state.mesh.dims, f64);
  //timer_heff=timer::start();

  for(unsigned i=0;i<Fieldterms.size();++i){
    solution+=Fieldterms[i]->h(state);
  }
  //time_heff+=af::timer::stop(timer_heff);
  return solution;
}

af::array NewLlg::fdmdt(const State& state){
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
  for(unsigned i=0;i<Fieldterms.size();++i){
    solution+=Fieldterms[i]->E(state);
  }
  return solution;
}

