#include "adaptive_runge_kutta.hpp"

AdaptiveRungeKutta::AdaptiveRungeKutta(callback_function f_in, std::string scheme_in, Controller controller): 
  f(f_in), scheme (scheme_in), controller(controller)
    {
    if (scheme == "RKF45") {
        std::cout << "Integrators: Initializing RKF45 method." << std::endl;
    }
    else if (scheme == "DP45") {
        std::cout << "Integrators: Initializing DP45 method." << std::endl;
    }
    else if (scheme == "BS45") {
        std::cout << "Integrators: Initializing BS45 method." << std::endl;
    }
    else {
        std::cout<< "Error: Integration method not found, please check the documantation" << std::endl;
        exit (EXIT_FAILURE);
    }
}

void AdaptiveRungeKutta::step(af::array& m, double& t){
    af::array mtemp;
    do{
        if (scheme == "RKF45") {
            mtemp = RKF45(m, t, h, err);
        }
        else if (scheme == "DP45")  {
            mtemp = DP45(m, t, h, err);
        }
        else {
            mtemp = BS45(m, t, h, err); 
        }
    }while(!controller.success(err,h));

    t+=h; //h is the actual timestep taken by the controller
    h=controller.get_hnext();
    m+=mtemp;
}

// Runge-Kutta-Fehlberg method with stepsize control
af::array AdaptiveRungeKutta::RKF45(const af::array& m, const double t, const double dt, double& err)
{
    af::array k1   =  dt * f(t             , m                                                                                                    );
    af::array k2   =  dt * f(t + 1./4.*dt  , m   +    1./4.    * k1                                                                               );
    af::array k3   =  dt * f(t + 3./8.*dt  , m   +    3./32.   * k1  + 9/32.       * k2                                                           );
    af::array k4   =  dt * f(t + 12./13.*dt, m   + 1932./2197. * k1  - 7200./2197. * k2   +  7296./2197. * k3                                     );
    af::array k5   =  dt * f(t + dt        , m   +  439./216.  * k1  -     8.      * k2   +  3680./513.  * k3  -   845./4104. * k4                );
    af::array k6   =  dt * f(t + 1./2.*dt  , m   -    8./27.   * k1  +     2.      * k2   -  3544./2565. * k3  +  1859./4104. * k4  - 11./40. * k5);
  
    af::array sumbk = 16./135. * k1 + 6656./12825.* k3 + 28561./56430.* k4 -9./50. * k5 + 2./55. *k6;
    af::array rk_error = sumbk - ( 25./216. * k1 + 1408./2565. * k3 + 2197./4104. * k4 -1./5. * k5);
  
    //rk_abs_error = maxnorm(rk_error);
    err=maxnorm(rk_error/controller.givescale(max(m,m+sumbk)));
    return sumbk;
}

//TODO
// FOR DP and BS, check why error is rising at the beginning of analytical example and then decreases again, maybe use different starting values
// Dormand-Prince 4/5 method
af::array AdaptiveRungeKutta::DP45(const af::array& m, const double t, const double dt, double& err)
{
    // Iterating over a-matrix and calculating k[i]s
    af::array k[8];
    const int s = 7;

    double a[8][7]={{0}};
    double e[8]={0};
    double c[8]={0};

    c[2]=0.2,c[3]=0.3,c[4]=0.8,c[5]=8.0/9.0, c[6]=1, c[7]=1;
    e[1]=71.0/57600.0,e[3]=-71.0/16695.0,e[4]=71.0/1920.0,e[5]=-17253.0/339200.0,e[6]=22.0/525.0,e[7]=-1.0/40.0,
    a[2][1]=0.2,
    a[3][1]=3.0/40.0,a[3][2]=9.0/40.0,
    a[4][1]=44.0/45.0,a[4][2]=-56.0/15.0,a[4][3]=32.0/9.0,
    a[5][1]=19372.0/6561.0,a[5][2]=-25360.0/2187.0,a[5][3]=64448.0/6561.0,a[5][4]=-212.0/729.0,
    a[6][1]=9017.0/3168.0,a[6][2]=-355.0/33.0,a[6][3]=46732.0/5247.0,a[6][4]=49.0/176.0,a[6][5]=-5103.0/18656.0,
    a[7][1]=35.0/384.0,a[7][3]=500.0/1113.0,a[7][4]=125.0/192.0,a[7][5]=-2187.0/6784.0,a[7][6]=11.0/84.0;

    
    const bool llg_wasnormalized=true;//TODO
    if(controller.get_reject() || (( controller.get_counter_reject() + controller.get_counter_accepted()) <=0) || llg_wasnormalized){
        k[1]   =  f(t, m);
        }
    else{
        k[1]=k[s];
    }
    for(int i=2;i<=s;i++){
        af::array rktemp=af::constant(0.0, m.dims(0), m.dims(1), m.dims(2), m.dims(3), f64);
        for(int j=1;j<i;j++){
            rktemp+=a[i][j] * k[j];
        }
        rktemp*=dt;
        k[i]= f(t+c[i], m + rktemp); 
    }
    //Local extrapolation using 5th order approx
    af::array sumbk=af::constant(0.0, m.dims(0), m.dims(1), m.dims(2), m.dims(3), f64);
    for(int i=1;i<s;i++){
        sumbk+=a[s][i]*k[i];
    }
    sumbk*=dt;
    //Error estimation using 4th order approx
    af::array rk_error=af::constant(0.0, m.dims(0), m.dims(1), m.dims(2), m.dims(3), f64);
    for(int i=1;i<=s;i++){
        rk_error+=e[i]*k[i];
    }
    rk_error*=dt;
    //!!!Note: here e is already the difference between the ususal b and bhat!!!! (no rk_error=sumbk-rk_error)
    err=maxnorm(rk_error/controller.givescale(max(m,m+sumbk)));
    return sumbk;
}


// Bogacki 4,5 method with sigle error andstepsize control
af::array AdaptiveRungeKutta::BS45(const af::array& m, const double t, const double dt , double& err)
{
    af::array k[9];
    const int s=8;
    double a[9][8]={{0.}};
    double b[9]={0.};
    double c[9]={0.};

    a[2][1] = 1.0e0/6.0e0;
    a[3][1] = 2.e0/27.e0;
    a[3][2] = 4.e0/27.e0;
    a[4][1] = 183.e0/1372.e0;
    a[4][2] = -162.e0/343.e0;
    a[4][3] = 1053.e0/1372.e0;
    a[5][1] = 68.e0/297.e0;
    a[5][2] = -4.e0/11.e0;
    a[5][3] = 42.e0/143.e0;
    a[5][4] = 1960.e0/3861.e0;
    a[6][1] = 597.e0/22528.e0;
    a[6][2] = 81.e0/352.e0;
    a[6][3] = 63099.e0/585728.e0;
    a[6][4] = 58653.e0/366080.e0;
    a[6][5] = 4617.e0/20480.e0;
    a[7][1] = 174197.e0/959244.e0;
    a[7][2] = -30942.e0/79937.e0;
    a[7][3] = 8152137.e0/19744439.e0;
    a[7][4] = 666106.e0/1039181.e0;
    a[7][5] = -29421.e0/29068.e0;
    a[7][6] = 482048.e0/414219.e0;
    a[8][1] = 587.e0/8064.e0;
    a[8][2] = 0.e0;
    a[8][3] = 4440339.e0/15491840.e0;
    a[8][4] = 24353.e0/124800.e0;
    a[8][5] = 387.e0/44800.e0;
    a[8][6] = 2152.e0/5985.e0;
    a[8][7] = 7267.e0/94080.e0;
  //C
  //C
    c[1] = 0.e0;
    c[2] = 1.e0/6.e0;
    c[3] = 2.e0/9.e0;
    c[4] = 3.e0/7.e0;
    c[5] = 2.e0/3.e0;
    c[6] = 3.e0/4.e0;
    c[7] = 1.e0;
    c[8] = 1.e0;
  //C  The coefficients B(*) refer to the formula of order 4.
    b[1] = 2479.e0/34992.e0;
    b[2] = 0.e0;
    b[3] = 123.e0/416.e0;
    b[4] = 612941.e0/3411720.e0;
    b[5] = 43.e0/1440.e0;
    b[6] = 2272.e0/6561.e0;
    b[7] = 79937.e0/1113912.e0;
    b[8] = 3293.e0/556956.e0;
    const bool llg_wasnormalized=true;//TODO
    if(controller.get_reject() || (( controller.get_counter_reject() + controller.get_counter_accepted()) <=0) || llg_wasnormalized){
    //if(reject || calls==0 || llg_wasnormalized){
    //Note: in generalized rkcall use: if(FSAL==false || reject || calls==0 || llg_wasnormalized){
      k[1]   =  f(t, m);
      }
    else
      k[1]=k[s];
    for(int i=2;i<=s;i++){
        af::array rktemp=af::constant(0.0, m.dims(0), m.dims(1), m.dims(2), m.dims(3), f64);
      for(int j=1;j<i;j++){
        rktemp+=a[i][j] * k[j];
      }
      rktemp*=dt;
      k[i]= f(t + c[i], m + rktemp); 
    }
  
    af::array sumbk=af::constant(0.0, m.dims(0), m.dims(1), m.dims(2), m.dims(3), f64);
    for(int i=1;i<s;i++){
      sumbk+=a[s][i]*k[i];
    }
    sumbk*=dt;
  
  
    af::array rk_error=af::constant(0.0, m.dims(0), m.dims(1), m.dims(2), m.dims(3), f64);
    for(int i=1;i<=s;i++){
      rk_error+=b[i]*k[i];
    }
    rk_error*=dt;
    rk_error=sumbk-rk_error;
    err=maxnorm(rk_error/controller.givescale(max(m,m+sumbk)));
  
    return sumbk;
}
