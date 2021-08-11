#pragma once
#include "arrayfire.h"

namespace magnumaf {

namespace ads_control {
inline af::array givescale(const af::array& val, double atol, double rtol) { return atol + rtol * af::abs(val); }
inline double givescale(double val, double atol, double rtol) { return atol + rtol * std::abs(val); }
} // namespace ads_control

class Controller {
  public:
    bool success(const double err,
                 double& h); // Decide whether step is acceppted

    Controller(double hmin = 1e-15, double hmax = 3.5e-10, double atol = 1e-6, double rtol = 1e-6, bool verbose = false)
        : hmin_(hmin), hmax_(hmax), atol_(atol), rtol_(rtol), verbose_(verbose) {}

    double get_hnext() const { return hnext_; }; // # of rejections
    bool get_reject() const { return reject_; };
    af::array givescale(const af::array& a) const { return ads_control::givescale(a, atol_, rtol_); };
    af::array givescale(double a) const { return ads_control::givescale(a, atol_, rtol_); };
    // af::array givescale(const af::array& a) const { return atol + rtol * af::abs(a); };

    // Access counters in read only
    unsigned long long int counter_reject() const { return counter_reject_; };     // # of rejections
    unsigned long long int counter_accepted() const { return counter_accepted_; }; // # of accepced steps
    unsigned long long int counter_hmax() const { return counter_hmax_; };         // # of rejections
    unsigned long long int counter_hmin() const { return counter_hmin_; };         // # of rejections
    unsigned long long int counter_maxscale() const { return counter_maxscale_; }; // # of rejections
    unsigned long long int counter_minscale() const { return counter_minscale_; }; // # of rejections

  private:
    const double hmin_;
    const double hmax_;
    // Scale function return= atol + abs(y) * rtol
    const double atol_;   // Tolerated absolute error
    const double rtol_;   // Tolerated relative error
    bool verbose_{false}; // Switch verbose mode

    // Numerical Recipies 3rd Edition suggests these values:
    const double beta_ = 0.4 / 5.0;
    const double alpha_ = 0.2 - 0.75 * beta_;
    const double headroom_ = 0.9;
    const double minscale_ = 0.2;
    const double maxscale_ = 10.;

    bool reject_{false};
    double errold_{1.0e-4}; // This value is max(err, 1.0e-4)
    double hnext_{0.0};

    // Counters to check performance
    unsigned long long int counter_reject_{0};   // # of rejections
    unsigned long long int counter_accepted_{0}; // # of accepced steps
    unsigned long long int counter_hmax_{0};     // # of rejections
    unsigned long long int counter_hmin_{0};     // # of rejections
    unsigned long long int counter_maxscale_{0}; // # of rejections
    unsigned long long int counter_minscale_{0}; // # of rejections
    // Member of LLG:
    // double  err{0};      // Estimated error
};
} // namespace magnumaf
