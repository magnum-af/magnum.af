#include "lbfgs_minimizer.hpp"

LBFGS_Minimizer::LBFGS_Minimizer()
{
}

// Calculation of effective field
af::array LBFGS_Minimizer::_h(const State& state){
    if( _llgterms.size() == 0){
        std::cout<<"ERROR: minimizer.cpp: Number of _llgterms == 0. Please add at least one term to LBFGS_Minimizer._llgterms! Aborting..."<<std::endl;
        exit (EXIT_FAILURE);
     }
    af::timer timer=af::timer::start();
    af::array solution = _llgterms[0]->h(state);
    for(unsigned i = 1; i < _llgterms.size() ; ++i ){
        solution+=_llgterms[i]->h(state);
    }
    _time_h+=af::timer::stop(timer);
    return solution;
}

void LBFGS_Minimizer::minimize(State& state){
    af::timer timer = af::timer::start();
    af::print("h in minimize", af::mean(_h(state)));//TODEL
    std::cout << "LBFGS_Minimizer: minimize in [s]: " << af::timer::stop(timer) << std::endl;
}; 
