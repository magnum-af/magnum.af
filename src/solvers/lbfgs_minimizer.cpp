#include "lbfgs_minimizer.hpp"

LBFGS_Minimizer::LBFGS_Minimizer()
{
}

// Calculation of effective field
af::array LBFGS_Minimizer::h(const State& state){
    if( llgterms.size() == 0){
        std::cout<<"ERROR: minimizer.cpp: Number of llgterms == 0. Please add at least one term to LBFGS_Minimizer.llgterms! Aborting..."<<std::endl;
        exit (EXIT_FAILURE);
     }
    af::timer timer=af::timer::start();
    af::array solution = llgterms[0]->h(state);
    for(unsigned i = 1; i < llgterms.size() ; ++i ){
        solution+=llgterms[i]->h(state);
    }
    time_h+=af::timer::stop(timer);
    return solution;
}

void LBFGS_Minimizer::minimize(State& state){
    af::timer timer = af::timer::start();
    af::print("h in minimize", af::mean(h(state)));//TODO

    std::cout << "LBFGS_Minimizer: time = " << af::timer::stop(timer) << std::endl;
}; 
