#ifndef CONTROLLER_H
#define CONTROLLER_H
#include "arrayfire.h"
#include <algorithm>
#include <iostream>


class Controller{
    public:
        bool success(const double err, double& h);//Decide whether step is acceppted

        const double hmin{1e-15};
        const double hmax{3.5e-10};
        double get_hnext () const { return hnext ;};// # of rejections
        af::array  givescale(const af::array& a){return atol+rtol*af::abs(a);}; 
        //Access counters in read only 
        unsigned int get_counter_reject  () const { return counter_reject  ;};// # of rejections
        unsigned int get_counter_accepted() const { return counter_accepted;};// # of accepced steps
        unsigned int get_counter_hmax    () const { return counter_hmax    ;};// # of rejections
        unsigned int get_counter_hmin    () const { return counter_hmin    ;};// # of rejections
        unsigned int get_counter_maxscale() const { return counter_maxscale;};// # of rejections
        unsigned int get_counter_minscale() const { return counter_minscale;};// # of rejections
    private:
        //Numerical Recipies 3rd Edition suggests these values:
        const double beta = 0.4/5.0;
        const double alpha = 0.2 - 0.75*beta;
        const double headroom = 0.9;
        const double minscale = 0.2;
        const double maxscale = 10.;

        // Scale function return= atol + abs(y) * rtol
        const double atol{1e-8};//Tolerated absolute error
        const double rtol{1e-8};//Tolerated relative error

        bool reject{false};
        double errold{1.0e-4};//This value is max(err,1.0e-4)
        double hnext;

        //Counters to check performance
        unsigned int counter_reject{0};// # of rejections
        unsigned int counter_accepted{0};// # of accepced steps
        unsigned int counter_hmax{0};// # of rejections
        unsigned int counter_hmin{0};// # of rejections
        unsigned int counter_maxscale{0};// # of rejections
        unsigned int counter_minscale{0};// # of rejections
    //Member of LLG:
    //double  err{0};      // Estimated error 

};


#endif
