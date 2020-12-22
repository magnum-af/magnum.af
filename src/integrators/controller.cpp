#include "controller.hpp"
#include "util/color_string.hpp" // for red
#include <cmath>

namespace magnumafcpp {

bool Controller::success(const double err, double& h) {
    double scale;

    if (err <= 1.0) {
        if (err == 0.0) {
            scale = maxscale_;
        } else {
            scale = headroom_ * std::pow(err, -alpha_) * std::pow(errold_, beta_);
            if (scale < minscale_) {
                scale = minscale_;
                counter_maxscale_++;
            }
            if (scale > maxscale_) {
                scale = maxscale_;
                counter_minscale_++;
            }
        }
        if (reject_) {
            hnext_ = h * std::min(scale, 1.0); // Do not increase stepsize if
                                               // previous try was rejected
        } else {
            hnext_ = h * scale;
        }
        if (hnext_ <= hmin_) {
            hnext_ = hmin_;
            counter_hmin_++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: proposed "
                   "stepsize hnext=%e <= minimum allowed stepsize "
                   "hmin=%e\nStepsize hmin will be forced and given error "
                   "bounds may become invalid.\n",
                   color_string::red("Warning:").c_str(), hnext_, hmin_);
        }
        if (hnext_ >= hmax_) {
            hnext_ = hmax_;
            counter_hmax_++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: proposed "
                   "stepsize hnext=%e >= maximum allowed stepsize "
                   "hmax=%e\nStepsize will be limited to hmax to preserve "
                   "given error bounds.\n",
                   color_string::red("Warning:").c_str(), hnext_, hmax_);
        }
        errold_ = std::max(err, 1.0e-4); // Why?
        reject_ = false;
        counter_accepted_++;
        return true;
    } else {
        scale = std::max(headroom_ * std::pow(err, -alpha_), minscale_);
        h *= scale;
        if (h <= hmin_) {
            h = hmin_;
            counter_hmin_++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: hmin=%e "
                   "reached in else, error bounds may be invalid.\n",
                   color_string::red("Warning:").c_str(), hmin_);
            return true;
        }
        if (h >= hmax_) {
            h = hmax_;
            counter_hmax_++;
            printf("%s Runge Kutta Adaptive Stepsize Controller: hmax=%e "
                   "reached in else.\n Stepsize is limited to hmax to preserve "
                   "error bounds.\n",
                   color_string::red("Warning:").c_str(), hmax_);
            return true;
        }
        reject_ = true;
        counter_reject_++;
        return false;
    }
}
} // namespace magnumafcpp
