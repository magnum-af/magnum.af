#include "controller.hpp"

bool Controller::success(const double err, double& h){
    double scale;
  
    if (err <= 1.0) {
        if (err == 0.0) {
            scale = maxscale;
        }
        else {
            scale=headroom*pow(err,-alpha)*pow(errold,beta);
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
            hnext=h*std::min(scale,1.0);//Do not increase stepsize if previous try was rejected
        }
        else {
            hnext=h*scale;
        }
        if (hnext<=hmin) {
            hnext=hmin;
            counter_hmin++;
            std::cout << red("Warning: hnext reached hmin in if, error bounds may be invalid")<<std::endl;
        }
        if (hnext>=hmax) {
            hnext=hmax;
            counter_hmax++;
            std::cout << red("Warning: hnext reached hmax in if, error bounds may be invalid") <<std::endl;
        }
        errold=std::max(err,1.0e-4);//Why?
        reject=false;
        counter_accepted++;
        return true;
    }
    else {
        scale=std::max(headroom*pow(err,-alpha),minscale);
        h *= scale;
        if (h<=hmin) {
            h=hmin;
            counter_hmin++;
            std::cout << red("Warning: hmin reached in else, error bounds may be invalid")<<std::endl;
            return true;
        }
        if (h>=hmax) {
            h=hmax;
            counter_hmax++;
            std::cout << red("Warning: hmax reached in else, error bounds may be invalid")<<std::endl;
            return true;
        }
        reject=true;
        counter_reject++;
        return false;
    }
}
