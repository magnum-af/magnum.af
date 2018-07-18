#include "adaptive_runge_kutta.hpp"

//af::array givem (const af::array& m){ return m;}//TODO del

AdaptiveRungeKutta::AdaptiveRungeKutta(callback_function f_in, std::string scheme_in, Controller controller, bool is_renormalize): 
  f(f_in), scheme (scheme_in), controller(controller), is_renormalize(is_renormalize)
    {
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

af::array AdaptiveRungeKutta::step(const af::array& m, const double t, double& dt){
    af::array mtemp;
    do{
        if (scheme == "RKF45") {
            mtemp = RKF45(m, t, h, err);//TODO
        }
        else { //if (scheme == "DP45") 
        }
    }while(!controller.success(err,h));

    dt=h; //dt is the actual timestep taken by the controller
    h=controller.get_hnext();
    mtemp+=m;

    //TODO if (state.Ms.isempty()) return  renormalize(mtemp);
    //return (renormalize_handle_zero_values(mtemp));//Normalized array where all initial values == 0 are set to 0
    if(is_renormalize) return renormalize(mtemp);
    else return mtemp;//TODO renorm
}

        //af::array (*foo)(const af::array&);
        //foo = &givem;

// Runge-Kutta-Fehlberg method with stepsize control
af::array AdaptiveRungeKutta::RKF45(const af::array& m, const double t, const double dt, double& err)
{
    af::array k1   =  dt * f(m                                                                                                    , t);
    af::array k2   =  dt * f(m   +    1./4.    * k1                                                                               , t + 1./4.*dt);
    af::array k3   =  dt * f(m   +    3./32.   * k1  + 9/32.       * k2                                                           , t + 3./8.*dt);
    af::array k4   =  dt * f(m   + 1932./2197. * k1  - 7200./2197. * k2   +  7296./2197. * k3                                     , t + 12./13.*dt);
    af::array k5   =  dt * f(m   +  439./216.  * k1  -     8.      * k2   +  3680./513.  * k3  -   845./4104. * k4                , t + dt);
    af::array k6   =  dt * f(m   -    8./27.   * k1  +     2.      * k2   -  3544./2565. * k3  +  1859./4104. * k4  - 11./40. * k5, t + 1./2.*dt);
  
    af::array sumbk     =                16./135.  * k1  +                       6656./12825.* k3  + 28561./56430.* k4    -9./50. * k5 + 2./55. *k6;
    af::array rk_error = sumbk - (      25./216.  * k1  +                       1408./2565. * k3  +  2197./4104. * k4    -1./5.  * k5             );
  
    //rk_abs_error = maxnorm(rk_error);
    err=maxnorm(rk_error/controller.givescale(max(m,m+sumbk)));
    return sumbk;
}

