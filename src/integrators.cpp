#include "integrators.hpp"

Integrator::Integrator(std::string scheme){

    if (scheme == "RKF45") {
        std::cout << "Integrators: Initializing RKF45 method." << std::endl;
        Controller controller=Controller();
    }

    if (scheme == "DP45") {
        std::cout << "Integrators: Initializing DP45 method." << std::endl;
        Controller controller=Controller();
    }
}

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
