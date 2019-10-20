#pragma once
#include "arrayfire.h"

namespace magnumaf{


class Controller{
    public:
        bool success(const float err, float& h);//Decide whether step is acceppted

        Controller(float hmin = 1e-15, float hmax = 3.5e-10, float atol = 1e-6, float rtol = 1e-6): hmin(hmin), hmax(hmax), atol(atol), rtol(rtol){}

        const float hmin;
        const float hmax;
        float get_hnext () const { return hnext ;};// # of rejections
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
        const float atol;//Tolerated absolute error
        const float rtol;//Tolerated relative error
    private:
        //Numerical Recipies 3rd Edition suggests these values:
        const float beta = 0.4/5.0;
        const float alpha = 0.2 - 0.75*beta;
        const float headroom = 0.9;
        const float minscale = 0.2;
        const float maxscale = 10.;


        bool reject{false};
        float errold{1.0e-4};//This value is max(err, 1.0e-4)
        float hnext;

        //Counters to check performance
        unsigned long long int counter_reject{0};// # of rejections
        unsigned long long int counter_accepted{0};// # of accepced steps
        unsigned long long int counter_hmax{0};// # of rejections
        unsigned long long int counter_hmin{0};// # of rejections
        unsigned long long int counter_maxscale{0};// # of rejections
        unsigned long long int counter_minscale{0};// # of rejections
    //Member of LLG:
    //float  err{0};      // Estimated error

};


}// namespace magnumaf
