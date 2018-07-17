#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include "arrayfire.h"
#include "controller.hpp"
#include "../func.hpp"

class Integrator {
    public: 
        Integrator(std::string scheme, bool is_set_call_time = false);
        af::array step(const af::array& m, const double dt);
        double get_call_time();
    private:
        const std::string scheme; // Denotes explicit integrator scheme s.a. RKF45 etc
        const bool is_set_call_time;
        double call_time; // Cummulative integrator call time
        Controller controller;
        af::array RKF45(af::array (*f)(af::array), const af::array& m, const double dt, double& err);
        //const int order;
};

//class RKF45 : public Integrator {
//    public: 
//        af::array step(const af::array& m, const double dt){
//            return af::constant(0,1,f64);//TODO remove 
//        }
//
//};

#endif



//#include "controller.hpp"
//#include <vector>
//class Integrator{
//    public:
//    Integrator(std::string scheme);
//
//    bool llg_wasnormalized{true};//set true after normalization in llgstep, false if not normalized
//                                 // DP45 uses this to decide wether k1 and heff has to be calculated
//                                 // -> leave the value true if normalization is on default
//    private:
//    int s; //Number of stages in RK method
//    bool reject{false};//True if last step was rejected TODO check when it is set true/false
//    bool FSAL=false;
//
//    // Counting step calls
//    unsigned int step_calls{0};
//    //TEMP
//    af::array k[14];
//    double a[14][14]={{0.}};
//    double b[14]={0.};
//    double bhat[14]={0.};
//    double c[14]={0.};
//    double e[14]={0.};
//    //TEMP
//
//
//    //Controller controller; // Initialize only in Adaptive Integrators = Controller(); 
//};

