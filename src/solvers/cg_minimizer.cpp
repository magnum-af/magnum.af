#include "cg_minimizer.hpp"

CG_Minimizer::CG_Minimizer()
{
}

// Calculation of effective field
af::array CG_Minimizer::h(const State& state){
    if( llgterms.size() == 0){
        std::cout<<"ERROR: minimizer.cpp: Number of llgterms == 0. Please add at least one term to CG_Minimizer.llgterms! Aborting..."<<std::endl;
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

void CG_Minimizer::minimize(State& state){
    af::timer timer = af::timer::start();
    af::print("h in minimize", af::mean(h(state)));//TODO

    std::cout << "CG_Minimizer: time = " << af::timer::stop(timer) << std::endl;
}; 
