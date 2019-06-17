#include "cg_minimizer.hpp"

CG_Minimizer::CG_Minimizer()
{
}

// Calculation of effective field
af::array CG_Minimizer::Heff(const State& state){
    if( llgterms_.size() == 0){
        std::cout<< bold_red("ERROR: LBFGS_Minimizer::Heff: Number of _llgterms == 0. Please add at least one term to LBFGS_Minimizer.llgterms_! Aborting...")<<std::endl;
        exit (EXIT_FAILURE);
     }
    af::timer timer=af::timer::start();
    af::array solution = llgterms_[0]->h(state);
    for(unsigned i = 1; i < llgterms_.size() ; ++i ){
        solution+=llgterms_[i]->h(state);
    }
    time_calc_heff_ += af::timer::stop(timer);
    return solution;
}

double CG_Minimizer::EnergyAndGradient(const State& state, af::array& gradient){
    if( llgterms_.size() == 0){
        std::cout<< bold_red("ERROR: LBFGS_Minimizer::Heff: Number of _llgterms == 0. Please add at least one term to LBFGS_Minimizer.llgterms_! Aborting...")<<std::endl;
        exit (EXIT_FAILURE);
     }
    af::timer timer=af::timer::start();
    //Avoiding array with zeros, starting loop with second term in llgterms
    af::array h = llgterms_[0]->h(state);
    double energy = llgterms_[0]->E(state, h);
    for(unsigned i = 1; i < llgterms_.size() ; ++i ){
        af::array temp_h = llgterms_[i]->h(state);
        h+=temp_h;
        energy+=llgterms_[i]->E(state,temp_h);
    }
    gradient = 1./(constants::mu0 * state.Ms) * cross4(state.m, cross4(state.m, h));
    time_calc_heff_ += af::timer::stop(timer);
    return energy;
}

void CG_Minimizer::Minimize(State& state){
    af::timer timer = af::timer::start();
    af::print("h in minimize", af::mean(Heff(state)));//TODO

    std::cout << "CG_Minimizer: time = " << af::timer::stop(timer) << std::endl;
};
