#include "adaptive_runge_kutta.hpp"

//af::array givem (const af::array& m){ return m;}//TODO del

AdaptiveRungeKutta::AdaptiveRungeKutta(callback_function f_in, std::string scheme_in): f(f_in), scheme (scheme_in){
    if (scheme == "RKF45") {
        std::cout << "Integrators: Initializing RKF45 method." << std::endl;
    }
    else if (scheme == "DP45") {
        std::cout << "Integrators: Initializing DP45 method." << std::endl;
    }
    else {
        std::cout<< "Error: Integration method not found, please check the documantation" << std::endl;
        exit (EXIT_FAILURE);
    }
}

af::array AdaptiveRungeKutta::step(const af::array& m, const double dt){
    if (scheme == "RKF45") {
        std::cout << "RKF45 step" << std::endl;
        //af::array (*foo)(const af::array&);
        //foo = &givem;
        double dummy;
        return RKF45(m, dt, dummy);//TODO
    }
    else { //if (scheme == "DP45") 
        return m;//TODO
    }
}


// Runge-Kutta-Fehlberg method with stepsize control
af::array AdaptiveRungeKutta::RKF45(const af::array& m, const double dt, double& err)
{
    af::array k1   =  dt * f(m                                                                                                    );
    af::array k2   =  dt * f(m   +    1./4.    * k1                                                                               );
    af::array k3   =  dt * f(m   +    3./32.   * k1  + 9/32.       * k2                                                           );
    af::array k4   =  dt * f(m   + 1932./2197. * k1  - 7200./2197. * k2   +  7296./2197. * k3                                     );
    af::array k5   =  dt * f(m   +  439./216.  * k1  -     8.      * k2   +  3680./513.  * k3  -   845./4104. * k4                );
    af::array k6   =  dt * f(m   -    8./27.   * k1  +     2.      * k2   -  3544./2565. * k3  +  1859./4104. * k4  - 11./40. * k5);
  
    af::array sumbk     =                16./135.  * k1  +                       6656./12825.* k3  + 28561./56430.* k4    -9./50. * k5 + 2./55. *k6;
    af::array rk_error = sumbk - (      25./216.  * k1  +                       1408./2565. * k3  +  2197./4104. * k4    -1./5.  * k5             );
  
    //rk_abs_error = maxnorm(rk_error);
    err=maxnorm(rk_error/controller.givescale(max(m,m+sumbk)));
    return sumbk;
}

