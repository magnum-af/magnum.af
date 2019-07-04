#pragma once
#include "arrayfire.h"

namespace magnumaf{


class Controller{
    public:
        bool success(const double err, double& h);//Decide whether step is acceppted

        Controller(double hmin = 1e-15, double hmax = 3.5e-10, double atol = 1e-6, double rtol = 1e-6): hmin(hmin), hmax(hmax), atol(atol), rtol(rtol){}

        const double hmin;
        const double hmax;
        double get_hnext () const { return hnext ;};// # of rejections
        af::array  givescale(const af::array& a){return atol+rtol*af::abs(a);};
        //Access counters in read only
        unsigned long long int get_counter_reject  () const { return counter_reject  ;};// # of rejections
        unsigned long long int get_counter_accepted() const { return counter_accepted;};// # of accepced steps
        unsigned long long int get_counter_hmax    () const { return counter_hmax    ;};// # of rejections
        unsigned long long int get_counter_hmin    () const { return counter_hmin    ;};// # of rejections
        unsigned long long int get_counter_maxscale() const { return counter_maxscale;};// # of rejections
        unsigned long long int get_counter_minscale() const { return counter_minscale;};// # of rejections

        bool get_reject() const { return reject;};

        //TODO better solution?
        // Scale function return= atol + abs(y) * rtol
        const double atol;//Tolerated absolute error
        const double rtol;//Tolerated relative error
    private:
        //Numerical Recipies 3rd Edition suggests these values:
        const double beta = 0.4/5.0;
        const double alpha = 0.2 - 0.75*beta;
        const double headroom = 0.9;
        const double minscale = 0.2;
        const double maxscale = 10.;


        bool reject{false};
        double errold{1.0e-4};//This value is max(err, 1.0e-4)
        double hnext;

        //Counters to check performance
        unsigned long long int counter_reject{0};// # of rejections
        unsigned long long int counter_accepted{0};// # of accepced steps
        unsigned long long int counter_hmax{0};// # of rejections
        unsigned long long int counter_hmin{0};// # of rejections
        unsigned long long int counter_maxscale{0};// # of rejections
        unsigned long long int counter_minscale{0};// # of rejections
    //Member of LLG:
    //double  err{0};      // Estimated error

};


}// namespace magnumaf
