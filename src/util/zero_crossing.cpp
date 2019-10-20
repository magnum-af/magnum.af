#include <iostream>
#include <math.h>
#include "zero_crossing.hpp"

namespace magnumaf{


using namespace magnumaf;

ZeroCrossing::ZeroCrossing(std::function<float (float)> f, float precision, int max_runs, float ix_min, float ix_max, int ix_n, int verbose) :
    f(f), precision(precision), max_runs(max_runs), ix_min(ix_min), ix_max(ix_max), ix_n(ix_n), verbose(verbose)
{
}


std::array<float, 4> ZeroCrossing::run_loop(){
    float x_max_minus_sign = -1e300;
    float f_max_minus_sign = -1e300;

    float x_min_plus_sign = +1e300;
    float f_min_plus_sign = +1e300;

    float fx_prev = 0;
    for (int i = 0; i < ix_n + 1; i++){
        const float x = ix_min + (float)i/(float)ix_n * (ix_max - ix_min);
        const float fx = f(x);
        if(verbose > 2) std::cout << "x=" << x << ", f(x) = " << fx << std::endl;
        if(fx <= 0 and fx > f_max_minus_sign){
            f_max_minus_sign = fx;
            x_max_minus_sign = x;
        }
        if(fx >= 0 and fx < f_min_plus_sign){
            f_min_plus_sign = fx;
            x_min_plus_sign = x;
        }
        if (i > 0 and fx_prev > 0 and fx <= 0) break;
        if (i > 0 and fx_prev < 0 and fx >= 0) break;
        fx_prev = fx;
    }
    return std::array<float, 4> {x_max_minus_sign, f_max_minus_sign, x_min_plus_sign, f_min_plus_sign};
}


std::pair<float, float> ZeroCrossing::calc_x_and_f(){
    std::array<float, 4> result = {0};
    for (int i = 0; i < max_runs; i++){
        result = run_loop();
        // Updating x range for next run
        // Rescale x range if no crossing was found
        if (result[0] == -1e300 or result[2] == 1e300){
            if(verbose) std::cout << "No crossing found in given interval. Setting larger x range by a factor of 10" << std::endl;
            // Setting ix_min < 0 < ix_max and scale by factor 10
            if (ix_min > ix_max){
                float temp = ix_min;
                ix_min = ix_max;
                ix_max = temp;
            }
            if(ix_min > 0) ix_min *= -1;
            if(ix_max < 0) ix_max *= -1;
            ix_min *= 10;
            ix_max *= 10;
        }
        // Otherwise use obtained nearest values
        else{
            ix_min = result[0];
            ix_max = result[2];
        }
        if(ix_min == ix_max){
            if(verbose) std::cout << "Warning: ZeroCrossing: ix_min == ix_max = " << ix_min << ", break!" << std::endl;
            break;
        }
        if(std::min(fabs(result[1]), fabs(result[3])) < precision){
            if(verbose) std::cout << "Info: ZeroCrossing: precision reached" << std::endl;
            break;
        }
        if(verbose > 1) std::cout << "i = " << i << ", ix_min = " << ix_min << ", ix_max = " << ix_max << std::endl << std::endl;
    }
    // only return x and f(x) closer to 0:
    if(fabs(result[1]) < fabs(result[3])) return std::pair<float, float> (result[0], result[1]);
    else return std::pair<float, float> (result[2], result[3]);
}
}// namespace magnumaf
