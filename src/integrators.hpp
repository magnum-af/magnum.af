#ifndef INTEGRATORS_H
#define INTEGRATORS_H
#include "arrayfire.h"
#include "controller.hpp"
#include <vector>

class Integrator{
    public:
    Integrator(std::string scheme);

    bool llg_wasnormalized{true};//set true after normalization in llgstep, false if not normalized
                                 // DP45 uses this to decide wether k1 and heff has to be calculated
                                 // -> leave the value true if normalization is on default
    private:
    int s; //Number of stages in RK method
    bool reject{false};//True if last step was rejected TODO check when it is set true/false
    bool FSAL=false;

    // Counting step calls
    unsigned int step_calls{0};
    //TEMP
    af::array k[14];
    double a[14][14]={{0.}};
    double b[14]={0.};
    double bhat[14]={0.};
    double c[14]={0.};
    double e[14]={0.};
    //TEMP


    //Controller controller; // Initialize only in Adaptive Integrators = Controller(); 
};

#endif
