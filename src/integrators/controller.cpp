#include "controller.hpp"
#include "../misc.hpp"
#include <algorithm>
#include <stdio.h>


bool Controller::success(const double err, double& h){
    double scale;

    if (err <= 1.0) {
        if (err == 0.0) {
            scale = maxscale;
        }
        else {
            scale=headroom*pow(err, -alpha)*pow(errold, beta);
            if (scale < minscale) {
                scale=minscale;
                counter_maxscale++;
            }
            if (scale > maxscale) {
                scale=maxscale;
                counter_minscale++;
            }
        }
        if (reject) {
            hnext=h*std::min(scale, 1.0);//Do not increase stepsize if previous try was rejected
        }
        else {
            hnext=h*scale;
        }
        if (hnext<=hmin) {
            hnext=hmin;
            counter_hmin++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: proposed stepsize hnext=%e <= minimum allowed stepsize hmin=%e\nStepsize hmin will be forced and given error bounds may become invalid.\n", red("Warning:").c_str(), hnext, hmin);
        }
        if (hnext>=hmax) {
            hnext=hmax;
            counter_hmax++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: proposed stepsize hnext=%e >= maximum allowed stepsize hmax=%e\nStepsize will be limited to hmax to preserve given error bounds.\n", red("Warning:").c_str(), hnext, hmax);
        }
        errold=std::max(err, 1.0e-4);//Why?
        reject=false;
        counter_accepted++;
        return true;
    }
    else {
        scale=std::max(headroom*pow(err, -alpha), minscale);
        h *= scale;
        if (h<=hmin) {
            h=hmin;
            counter_hmin++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: hmin=%e reached in else, error bounds may be invalid.\n", red("Warning:").c_str(), hmin);
            return true;
        }
        if (h>=hmax) {
            h=hmax;
            counter_hmax++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: hmax=%e reached in else.\n Stepsize is limited to hmax to preserve error bounds.\n", red("Warning:").c_str(), hmax);
            return true;
        }
        reject=true;
        counter_reject++;
        return false;
    }
}
