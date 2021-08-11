#include "controller.hpp"
#include <cmath>

#include "util/color_string.hpp" // for red
#include <cmath>

namespace magnumaf {

void warning_hnext_lt_hmin(double hnext, double hmin) {
    printf("%s Runge Kutta Adaptive Stepsize Controller: proposed "
           "stepsize hnext=%e <= minimum allowed stepsize "
           "hmin=%e. Stepsize hmin will be forced and given error "
           "bounds may become invalid.\n",
           color_string::red("Warning:").c_str(), hnext, hmin);
}

void warning_h_lt_hmin(double h, double hmin) {
    printf("%s Runge Kutta Adaptive Stepsize Controller: tried step h=%e <= hmin=%e."
           "Error bounds may become invalid.\n",
           color_string::red("Warning:").c_str(), h, hmin);
}

void info_hnext_gt_hmax(double hnext, double hmax) {
    printf("%s Runge Kutta Adaptive Stepsize Controller: proposed "
           "stepsize hnext=%e >= maximum allowed stepsize "
           "hmax=%e\nStepsize will be limited to hmax.  This means that large timesteps are taken.\n",
           color_string::green("Info:").c_str(), hnext, hmax);
}

void info_h_gt_hmax(double h, double hmax) {
    printf("%s Runge Kutta Adaptive Stepsize Controller: tried step h=%e >= hmax=%e "
           ".\n Stepsize is limited to hmax. This means that large timesteps are taken.\n",
           color_string::red("Info:").c_str(), h, hmax);
}

bool Controller::success(const double err, double& h) {
    double scale = std::numeric_limits<double>::quiet_NaN();

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
            warning_hnext_lt_hmin(hnext_, hmin_);
            hnext_ = hmin_;
            counter_hmin_++;
        } else if (hnext_ >= hmax_) {
            if (verbose_) {
                info_hnext_gt_hmax(hnext_, hmax_);
            }
            hnext_ = hmax_;
            counter_hmax_++;
        }
        errold_ = std::max(err, 1.0e-4); // Why?
        reject_ = false;
        counter_accepted_++;
        return true;
    } else {
        scale = std::max(headroom_ * std::pow(err, -alpha_), minscale_);
        h *= scale;
        if (h <= hmin_) {
            warning_h_lt_hmin(h, hmin_);
            h = hmin_;
            counter_hmin_++;
            return true;
        } else if (h >= hmax_) {
            if (verbose_) {
                info_h_gt_hmax(h, hmax_);
            }
            h = hmax_;
            counter_hmax_++;
            return true;
        } else {
            reject_ = true;
            counter_reject_++;
            return false;
        }
    }
}
} // namespace magnumaf
