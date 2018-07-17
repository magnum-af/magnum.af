#include "integrator.hpp"

Integrator::Integrator(std::string scheme_in, bool is_set_call_time_in): scheme(scheme_in), is_set_call_time(is_set_call_time_in){

    if (scheme == "RKF45") {
        std::cout << "Integrators: Initializing RKF45 method." << std::endl;
        Controller controller=Controller();
    }
    else if (scheme == "DP45") {
        std::cout << "Integrators: Initializing DP45 method." << std::endl;
    }
    else {
        std::cout<< "Error: Integration method not found, please check the documantation" << std::endl;
    }
}


af::array Integrator::step(const af::array& m, const double dt){
    af::timer timer_integ;
    if (is_set_call_time) timer_integ = af::timer::start();
    if (scheme == "RKF45") {
        std::cout << "RKF45 step" << std::endl;
        //af::array out = RKF45(f, m, dt, 
    }
    if (is_set_call_time) call_time += af::timer::stop(timer_integ);

    return m;//TODO
}

// Runge-Kutta-Fehlberg method with stepsize control
af::array Integrator::RKF45(af::array (*f)(af::array), const af::array& m, const double dt, double& err)
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

//// Runge-Kutta-Fehlberg 5th order method
//array LLG::RKF5(const array& m, const double dt)
//{
//
//  heff =       fheff(m                                                                                                                         );
//  k1   =  dt * fdmdt(m                                                                                                                   , heff);
//  heff =       fheff(m   +    1./4.    * k1                                                                                                    );
//  k2   =  dt * fdmdt(m   +    1./4.    * k1                                                                                              , heff);
//  heff =       fheff(m   +    3./32.   * k1  + 9/32.       * k2                                                                                );
//  k3   =  dt * fdmdt(m   +    3./32.   * k1  + 9/32.       * k2                                                                          , heff);
//  heff =       fheff(m   + 1932./2197. * k1  - 7200./2197. * k2   +  7296./2197. * k3                                                          );
//  k4   =  dt * fdmdt(m   + 1932./2197. * k1  - 7200./2197. * k2   +  7296./2197. * k3                                                    , heff);
//  heff =       fheff(m   +  439./216.  * k1  -     8.      * k2   +  3680./513.  * k3  -   845./4104. * k4                                     );
//  k5   =  dt * fdmdt(m   +  439./216.  * k1  -     8.      * k2   +  3680./513.  * k3  -   845./4104. * k4                               , heff);
//  heff =       fheff(m   -    8./27.   * k1  +     2.      * k2   -  3544./2565. * k3  +  1859./4104. * k4  - 11./40. * k5                     );
//  k6   =  dt * fdmdt(m   -    8./27.   * k1  +     2.      * k2   -  3544./2565. * k3  +  1859./4104. * k4  - 11./40. * k5               , heff);
//
//  sumbk     =                16./135.  * k1  +                       6656./12825.* k3  + 28561./56430.* k4    -9./50. * k5 + 2./55. *k6         ;
//  return sumbk;
//}



//af::array Integrator::step(const af::array& m, const double dt, double& err){
//  
//    // Iterating over a-matrix and calculating k[i]s
//    if(reject || calls==0 || llg_wasnormalized){
//        array heff =  fheff(m);
//        k[1]   =  fdmdt(m, heff);
//    }
//    else{
//        k[1]=k[s];
//    }
//    for(int i=2;i<=s;i++){
//          rktemp=constant(0.0,state0.mesh.n0, state0.mesh.n1,state0.mesh.n2,3,f64);
//        for(int j=1;j<i;j++){
//            rktemp+=a[i][j] * k[j];
//        }
//        rktemp*=dt;
//        heff= fheff(m + rktemp); 
//        k[i]= fdmdt(m + rktemp,heff); 
//    }
//    //Local extrapolation using 5th order approx
//    sumbk=constant(0.0,state0.mesh.n0, state0.mesh.n1,state0.mesh.n2,3,f64);
//    for(int i=1;i<s;i++){
//        sumbk+=a[s][i]*k[i];
//    }
//    sumbk*=dt;
//    //Error estimation using 4th order approx
//    rk_error=constant(0.0,state0.mesh.n0, state0.mesh.n1,state0.mesh.n2,3,f64);
//    for(int i=1;i<=s;i++){
//        rk_error+=e[i]*k[i];
//    }
//    rk_error*=dt;
//    //!!!Note: here e is already the difference between the ususal b and bhat!!!! (no rk_error=sumbk-rk_error)
//    err=maxnorm(rk_error/givescale(max(m,m+sumbk)));
//    return sumbk;
//}



//// Dormand-Prince 4/5 method
//array LLG::DP45(const array& m, const double dt, double& err)
//{
//  // Iterating over a-matrix and calculating k[i]s
//  if(reject || calls==0 || llg_wasnormalized){
//    heff =  fheff(m);
//    k[1]   =  fdmdt(m, heff);
//    }
//  else
//    k[1]=k[s];
//  for(int i=2;i<=s;i++){
//      rktemp=constant(0.0,state0.mesh.n0, state0.mesh.n1,state0.mesh.n2,3,f64);
//    for(int j=1;j<i;j++){
//      rktemp+=a[i][j] * k[j];
//    }
//    rktemp*=dt;
//    heff= fheff(m + rktemp); 
//    k[i]= fdmdt(m + rktemp,heff); 
//  }
//  //Local extrapolation using 5th order approx
//  sumbk=constant(0.0,state0.mesh.n0, state0.mesh.n1,state0.mesh.n2,3,f64);
//  for(int i=1;i<s;i++){
//    sumbk+=a[s][i]*k[i];
//  }
//  sumbk*=dt;
//  //Error estimation using 4th order approx
//  rk_error=constant(0.0,state0.mesh.n0, state0.mesh.n1,state0.mesh.n2,3,f64);
//  for(int i=1;i<=s;i++){
//    rk_error+=e[i]*k[i];
//  }
//  rk_error*=dt;
//  //!!!Note: here e is already the difference between the ususal b and bhat!!!! (no rk_error=sumbk-rk_error)
//  err=maxnorm(rk_error/givescale(max(m,m+sumbk)));
//  return sumbk;
//}
