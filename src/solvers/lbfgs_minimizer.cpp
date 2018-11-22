#include "lbfgs_minimizer.hpp"

LBFGS_Minimizer::LBFGS_Minimizer()
{
}

// Calculation of effective field
af::array LBFGS_Minimizer::Heff(const State& state){
    if( llgterms_.size() == 0){
        std::cout<<"ERROR: LBFGS_Minimizer::Heff: Number of _llgterms == 0. Please add at least one term to LBFGS_Minimizer.llgterms_! Aborting..."<<std::endl;
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

af::array LBFGS_Minimizer::Gradient(const State& state){
    return - state.param.alpha*state.param.gamma/(1.+pow(state.param.alpha,2)) * cross4(state.m, cross4(state.m, Heff(state)));
}

void LBFGS_Minimizer::Minimize(State& state){
    af::timer timer = af::timer::start();
    //af::print("h in minimize", af::mean(Heff(state)));//TODEL
    //af::print("Gradient", Gradient(state));//TODEL
    auto grad = Gradient(state);

    double eps  = 1e-12; //TODO CL_DBL_EPSILON;
    double eps2 = sqrt(eps);
    double epsr = pow(eps,0.9);
    double tolf = 1e-12; //TODO this->settings_.gradTol;
    double tolf2 = sqrt(tolf);
    double tolf3 = pow(tolf,0.3333333333333333333333333);
    //double f = objFunc.both(x0, grad);// objFunc.both calcs Heff and E for not calculating Heff double
    // TODO for starting, repalce linesearch with return value 1.0



    std::cout << "LBFGS_Minimizer: minimize in [s]: " << af::timer::stop(timer) << std::endl;
}; 
