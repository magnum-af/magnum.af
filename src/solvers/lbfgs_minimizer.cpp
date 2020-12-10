#include "lbfgs_minimizer.hpp"
#include "arrayfire.h"
#include "equations.hpp"
#include "field_terms/field_term.hpp"
#include "state.hpp"
#include "util/func.hpp"
#include "util/misc.hpp"
#include <algorithm>
#include <iomanip>
#include <list>
#include <memory>

namespace magnumafcpp {

LBFGS_Minimizer::LBFGS_Minimizer(double tolerance, size_t maxIter, int verbose)
    : tolerance_(tolerance), maxIter_(maxIter), verbose_(verbose) {}

LBFGS_Minimizer::LBFGS_Minimizer(vec_uptr_FieldTerm llgterms, double tolerance, size_t maxIter, int verbose)
    : fieldterms_(std::move(llgterms)), tolerance_(tolerance), maxIter_(maxIter), verbose_(verbose) {}

void abort_if_size_is_zero(std::size_t size) {
    if (size == 0) {
        std::cout << bold_red("ERROR: LBFGS_Minimizer::Heff: Number of "
                              "_llgterms == 0. Please add at least one term to "
                              "LBFGS_Minimizer.fieldterms_! Aborting...")
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}

///< Calculate gradient as m x (m x Heff), which is proportional to energy-dissipation term of the LLG
af::array Gradient(const State& state, const af::array& heff) {
    return cross4(state.m, cross4(state.m, heff)); // Note: Works equally well as current grad, is not Ms dependent
}

// Alternatives
// af::array Gradient(const State& state, const af::array& heff) {
//    return 1. / (constants::mu0 * state.Ms) * cross4(state.m, cross4(state.m, heff)); // NOTE: state.Ms must not be 0!
//    return -equations::LLG_damping(1, state.m, cross4(state.m, heff)); // Also works, but is stuck at step 15 in
//}

std::pair<double, af::array> EnergyAndGradient(const State& state, const vec_uptr_FieldTerm& fieldterms) {
    abort_if_size_is_zero(fieldterms.size());
    const auto [heff, energy] = fieldterm::Heff_in_Apm_and_E(fieldterms, state);
    return {energy, Gradient(state, heff)};
}

double mydot(const af::array& a, const af::array& b) { return full_inner_product(a, b); }
double mynorm(const af::array& a) { return sqrt(mydot(a, a)); }

double mxmxhMax(const State& state, const vec_uptr_FieldTerm& fieldterms) {
    return max_4d_abs(Gradient(state, fieldterm::Heff_in_Apm(fieldterms, state)));
}

/// LBFGS minimizer from Thomas Schrefl's bvec code
// TODO Currently Minimize() fails when called second time, i.e. when m already
// relaxed linesearch rate is 0 even with changed zeeman field
double LBFGS_Minimizer::Minimize(State& state) const {
    std::cout.precision(24);
    af::timer timer = af::timer::start();
    // af::print("h in minimize", af::mean(Heff(state)));//TODEL
    // af::print("Gradient", Gradient(state));//TODEL

    double eps = 2.22e-16;
    double eps2 = sqrt(eps);
    double epsr = pow(eps, 0.9);
    // double tolerance_ = 1e-6; //TODO find value of this->settings_.gradTol;
    double tolf2 = sqrt(tolerance_);
    double tolf3 = pow(tolerance_, 0.3333333333333333333333333);
    auto [f, grad] = EnergyAndGradient(state, fieldterms_);

    // double f = objFunc.both(x0, grad);// objFunc.both calcs Heff and E for
    // not calculating Heff double
    // NOTE: objFunc.both calcs gradient and energy E
    double gradNorm = max_4d_abs(grad);
    // af::print("grad", grad);
    if (verbose_ > 0) {
        std::cout << "f= " << f << std::endl;
    }
    if (gradNorm < (epsr * (1 + fabs(f)))) {
        return f;
    }
    const size_t m = 5;

    // array container
    af::array af_zero = (af::constant(0.0, mesh::dims_v(state.mesh), f64));
    std::array<af::array, m> sVector{{af_zero, af_zero, af_zero, af_zero, af_zero}};
    std::array<af::array, m> yVector{{af_zero, af_zero, af_zero, af_zero, af_zero}};

    // vector container
    // std::vector<af::array> sVector; //= af::constant(0.0, mesh::dims_v(state.mesh),
    // f64); std::vector<af::array> yVector; //= af::constant(0.0,
    // mesh::dims_v(state.mesh), f64); for (size_t i = 0; i < m; i++){
    //    sVector.push_back(af::constant(0.0, mesh::dims_v(state.mesh), f64));
    //    yVector.push_back(af::constant(0.0, mesh::dims_v(state.mesh), f64));
    //}

    std::array<double, m> alpha = {0., 0., 0., 0., 0.};
    af::array q = af::constant(0.0, mesh::dims_v(state.mesh), f64);
    af::array s = af::constant(0.0, mesh::dims_v(state.mesh), f64);
    af::array y = af::constant(0.0, mesh::dims_v(state.mesh), f64);
    double f_old = 0.0;
    af::array grad_old = af::constant(0.0, mesh::dims_v(state.mesh), f64);
    af::array x_old = state.m;

    size_t iter = 0, globIter = 0;
    double H0k = 1;
    do {
        // TODO TODEL
        //
        // for (size_t i = 0; i < m; i++){
        //    printf("size svec = %i \n", sVector.size());
        //    af::print("svec", af::mean(sVector[i], 0));
        //}
        // END TODO
        int cgSteps = 0;
        // this->minIterCount_++;
        f_old = f;
        x_old = state.m;
        // af::print("state.m", af::mean(state.m, 0));
        grad_old = grad;

        q = grad;
        const int k = std::min(m, iter);
        for (int i = k - 1; i >= 0; i--) {
            const double rho = 1.0 / mydot(sVector[i], yVector[i]);
            alpha[i] = rho * mydot(sVector[i], q);
            q -= alpha[i] * yVector[i];
        }
        q = H0k * q; // NOTE: cg step skipped and only used else
        for (int i = 0; i < k; i++) {
            const double rho = 1.0 / mydot(sVector[i], yVector[i]);
            const double beta = rho * mydot(yVector[i], q);
            q += sVector[i] * (alpha[i] - beta);
        }

        // line 291 in .fe: decent = -grad.dot(q)
        double phiPrime0 = -mydot(grad, q);
        if (phiPrime0 > 0) {
            q = grad;
            iter = 0;
            if (verbose_ > 2) {
                std::cout << "descent " << std::endl;
            }
            phiPrime0 = -mydot(grad, q);
        }

        // TODO objFunc.updateTol(gradNorm);
        // TODO check version Schrefl vs Flo
        if (-mydot(grad, q) > -1e-15) {
            gradNorm = max_4d_abs(grad);
            if (gradNorm < eps * (1. + fabs(f)) and verbose_ > 0) {
                std::cout << "Minimizer: Convergence reached (due to almost "
                             "zero gradient (|g|="
                          << gradNorm << " < " << eps * (1. + fabs(f)) << ")!" << std::endl;
            }
        }
        // const double rate = MyMoreThuente<T, decltype(objFunc),
        // 1>::linesearch(f, x_old, x0, grad, -q, objFunc, tolf);
        const double rate = linesearch(state, f, x_old, grad, -q, tolerance_);
        // TODO/old/todel//const double rate = this->cvsrch(state, x_old, x0, f,
        // grad, -q, tolf);
        if (rate == 0.0) {
            if (verbose_ > 0) {
                std::cout << red("Warning: LBFGS_Minimizer: linesearch returned rate "
                                 "== 0.0. This should not happen, elaborate! (maybe "
                                 "m is already relaxed?)")
                          << std::endl;
                // std::cout << bold_red("Error: LBFGS_Minimizer: linesearch
                // failed, rate == 0.0") << std::endl;
            }
            // TODO//exit(0);// why exit with 0 (= sucess)? maybe change to
            // return(0) or, (if necessary)//exit(EXIT_FAILURE);
        }
        double f1 = 1 + fabs(f);
        gradNorm = max_4d_abs(grad);
        if (gradNorm < (epsr * f1)) {
            return f;
        }
        s = state.m - x_old;
        if (verbose_ > 1) {
            // std::cout <<
            // std::setprecision(std::numeric_limits<double>::digits10 + 1);
            std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 12);
            std::cout << "bb> " << globIter << " " << f << " "
                      << " " << gradNorm << " " << cgSteps << " " << rate << std::endl;
        }

        if (of_convergence_.is_open()) {
            of_convergence_ << std::setprecision(std::numeric_limits<double>::digits10 + 12);
            of_convergence_ << (f_old - f) / (tolerance_ * f1) << "\t"
                            << max_4d_abs(s) / (tolf2 * (1 + max_4d_abs(state.m))) << "\t" << gradNorm / (tolf3 * f1)
                            << std::endl;
        }

        if (((f_old - f) < (tolerance_ * f1)) && (max_4d_abs(s) < (tolf2 * (1 + max_4d_abs(state.m)))) &&
            (gradNorm <= (tolf3 * f1))) {
            break;
        }
        y = grad - grad_old;
        double ys = mydot(y, s);
        if (ys <= eps2 * mynorm(y) * mynorm(s)) { // Dennis, Schnabel 9.4.1
            if (verbose_ > 2) {
                std::cout << iter << red("WARNING: LBFGS_Minimizer:: skipping update!") << std::endl;
            }
        } else {
            if (iter < m) {
                // af::print("s", af::mean(s, 0));
                sVector[iter] = s;
                // af::print("svec", af::mean(sVector[iter], 0));
                yVector[iter] = y;
            } else {
                std::rotate(sVector.begin(), sVector.begin() + 1, sVector.end());
                sVector[m - 1] = s;
                std::rotate(yVector.begin(), yVector.begin() + 1, yVector.end());
                yVector[m - 1] = y;
            }
            H0k = ys / mydot(y, y);
            iter++;
        }

        globIter++;

    } while (globIter < this->maxIter_);
    if (globIter >= this->maxIter_ && this->maxIter_ > 99 and verbose_ > 0) {
        std::cout << "WARNING : maximum number of iterations exceeded in LBFGS" << std::endl;
    }
    if (verbose_ > 0) {
        State temp = state;
        temp.m = grad;
        auto mxh = mxmxhMax(temp, fieldterms_);
        auto deltaF = f - f_old;
        std::cout << "aa> " << globIter << " " << deltaF << " "
                  << " " << mxh << std::endl;
    }
    return f;

    if (verbose_ > 0) {
        std::cout << "LBFGS_Minimizer: minimize in [s]: " << af::timer::stop(timer) << std::endl;
    }
}

double LBFGS_Minimizer::linesearch(State& state, double& fval, const af::array& x_old, af::array& g,
                                   const af::array& searchDir, const double tolf) const {
    double rate = 1.0;
    cvsrch(state, x_old, fval, g, rate, searchDir, tolf);
    return rate;
}

// static int cvsrch(P &objFunc, const vex::vector<double> &wa,
// vex::vector<double> &x, double &f, vex::vector<double> &g, double &stp, const
// vex::vector<double> &s, double tolf) {
// TODO//TODEL//int LBFGS_Minimizer::cvsrch(const State& state, const af::array
// &wa, af::array &x, double &f, af::array &g, const af::array &s, double tolf)
// {// ak = 1.0 == double &stp moved into function
int LBFGS_Minimizer::cvsrch(State& state, const af::array& wa, double& f, af::array& g, double& stp, const af::array& s,
                            const double tolf) const {
    // we rewrite this from MIN-LAPACK and some MATLAB code

    int info = 0;
    int infoc = 1;

    //  const int n        = x.dims(0)*x.dims(1)*x.dims(2)*x.dims(3); (void) n;
    const double xtol = 1e-15;
    const double ftol = 1.0e-4; // c1
    const double gtol = 0.9;    // c2
    const double eps = tolf;
    const double stpmin = 1e-15;
    const double stpmax = 1e15;
    const double xtrapf = 4;
    const int maxfev = 20;
    int nfev = 0;

    double dginit = mydot(g, s);
    if (dginit >= 0.0) {
        // no descent direction
        // TODO: handle this case
        std::cout << red("WARNING: LBFGS_Minimizer:: no descent ") << dginit << std::endl;
        return -1;
    }

    bool brackt = false;
    bool stage1 = true;

    double finit = f;
    double dgtest = ftol * dginit;
    double width = stpmax - stpmin;
    double width1 = 2 * width;
    // vex::vector<double> wa(x.queue_list(), x.size());
    // wa = x;

    double stx = 0.0;
    double fx = finit;
    double dgx = dginit;
    double sty = 0.0;
    double fy = finit;
    double dgy = dginit;

    double stmin;
    double stmax;

    while (true) {

        // make sure we stay in the interval when setting min/max-step-width
        if (brackt) {
            stmin = std::min(stx, sty);
            stmax = std::max(stx, sty);
        } else {
            stmin = stx;
            stmax = stp + xtrapf * (stp - stx);
        }

        // Force the step to be within the bounds stpmax and stpmin.
        stp = std::max(stp, stpmin);
        stp = std::min(stp, stpmax);

        // Oops, let us return the last reliable values
        if ((brackt && ((stp <= stmin) | (stp >= stmax))) | (nfev >= maxfev - 1) | (infoc == 0) |
            (brackt & (stmax - stmin <= xtol * stmax))) {
            stp = stx;
            if (verbose_ > 0)
                std::cout << "NOTE: LBFGS_Minimizer:: Oops, stp= " << stp << std::endl;
        }

        //// test new point
        //// x = wa + stp * s;
        //// objFunc.normalizeVector(x);
        // objFunc.update(stp, wa, s, x);
        state.m = wa + stp * s; // TODO check// this should be equivalent to
                                // objFunc.update(stp, wa, s, x);
        state.m = normalize_handle_zero_vectors(state.m);
        std::tie(f, g) = EnergyAndGradient(state, fieldterms_);
        // Note: using struc-binding via 'auto [f, g] = ' is a bug
        // as it would declare f,g in this while scope and shadow outer f,g
        nfev++;
        double dg = mydot(g, s);
        double ftest1 = finit + stp * dgtest;
        double ftest2 = finit + eps * fabs(finit);
        double ft = 2 * ftol - 1;

        // all possible convergence tests
        if ((brackt & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0))
            info = 6;

        // if ((stp == stpmax) & (f <= ftest1) & (dg <= dgtest))
        if ((stp == stpmax) & (f <= ftest2) & (dg <= dgtest))
            info = 5;

        // if ((stp == stpmin) & ((f > ftest1) | (dg >= dgtest)))
        if ((stp == stpmin) & ((f > ftest2) | (dg >= dgtest)))
            info = 4;

        if (nfev >= maxfev)
            info = 3;

        if (brackt & (stmax - stmin <= xtol * stmax))
            info = 2;

        if ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit)))
            info = 1;

        if (((f <= ftest2) && (ft * dginit >= dg)) && (fabs(dg) <= gtol * (-dginit)))
            info = 1;

        // terminate when convergence reached
        if (info != 0) {
            return -1;
        }

        // if (stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol)*dginit))
        if (stage1 & ((f <= ftest2) && (ft * dginit >= dg)) & (dg >= std::min(ftol, gtol) * dginit)) // approx wolfe 1
            stage1 = false;

        // if (stage1 & (f <= fx) & (f > ftest1)) {
        if (stage1 & (f <= fx) & (not((f <= ftest2) && (ft * dginit >= dg)))) { // not wolfe 1 --> not approx wolfe 1
            double fm = f - stp * dgtest;
            double fxm = fx - stx * dgtest;
            double fym = fy - sty * dgtest;
            double dgm = dg - dgtest;
            double dgxm = dgx - dgtest;
            double dgym = dgy - dgtest;

            int rsl = cstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc);
            (void)rsl;

            fx = fxm + stx * dgtest;
            fy = fym + sty * dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
        } else {
            // this is ugly and some variables should be moved to the class
            // scope
            int rsl = cstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc);
            (void)rsl;
        }

        if (brackt) {
            if (fabs(sty - stx) >= 0.66 * width1)
                stp = stx + 0.5 * (sty - stx);
            width1 = width;
            width = fabs(sty - stx);
        }
    }

    std::cout << bold_red("ERROR: LBFGS_Minimizer: cvsrch: XXXXXXXXXXXXXXXXXXXXXXXXXXXX "
                          "why I am here ? XXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxx ")
              << std::endl;
    exit(0);
    return 0;
}

int LBFGS_Minimizer::cstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, double& stp,
                           double& fp, double& dp, bool& brackt, double& stpmin, double& stpmax, int& info) const {
    info = 0;
    bool bound = false;

    // Check the input parameters for errors.
    if ((brackt & ((stp <= std::min(stx, sty)) | (stp >= std::max(stx, sty)))) | (dx * (stp - stx) >= 0.0) |
        (stpmax < stpmin)) {
        return -1;
    }

    double sgnd = dp * (dx / fabs(dx));

    double stpf = 0;
    double stpc = 0;
    double stpq = 0;

    if (fp > fx) {
        info = 1;
        bound = true;
        double theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
        double s = std::max(theta, std::max(dx, dp));
        double gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if (stp < stx)
            gamma = -gamma;
        double p = (gamma - dx) + theta;
        double q = ((gamma - dx) + gamma) + dp;
        double r = p / q;
        stpc = stx + r * (stp - stx);
        stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.) * (stp - stx);
        if (fabs(stpc - stx) < fabs(stpq - stx))
            stpf = stpc;
        else
            stpf = stpc + (stpq - stpc) / 2;
        brackt = true;
    } else if (sgnd < 0.0) {
        info = 2;
        bound = false;
        double theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        double s = std::max(theta, std::max(dx, dp));
        double gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
        if (stp > stx)
            gamma = -gamma;

        double p = (gamma - dp) + theta;
        double q = ((gamma - dp) + gamma) + dx;
        double r = p / q;
        stpc = stp + r * (stx - stp);
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (fabs(stpc - stp) > fabs(stpq - stp))
            stpf = stpc;
        else
            stpf = stpq;
        brackt = true;
    } else if (fabs(dp) < fabs(dx)) {
        info = 3;
        bound = 1;
        double theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        double s = std::max(theta, std::max(dx, dp));
        double gamma = s * sqrt(std::max(0., (theta / s) * (theta / s) - (dx / s) * (dp / s)));
        if (stp > stx)
            gamma = -gamma;
        double p = (gamma - dp) + theta;
        double q = (gamma + (dx - dp)) + gamma;
        double r = p / q;
        if ((r < 0.0) & (gamma != 0.0)) {
            stpc = stp + r * (stx - stp);
        } else if (stp > stx) {
            stpc = stpmax;
        } else {
            stpc = stpmin;
        }
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (brackt) {
            if (fabs(stp - stpc) < fabs(stp - stpq)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
        } else {
            if (fabs(stp - stpc) > fabs(stp - stpq)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
        }
    } else {
        info = 4;
        bound = false;
        if (brackt) {
            double theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
            double s = std::max(theta, std::max(dy, dp));
            double gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
            if (stp > sty)
                gamma = -gamma;

            double p = (gamma - dp) + theta;
            double q = ((gamma - dp) + gamma) + dy;
            double r = p / q;
            stpc = stp + r * (sty - stp);
            stpf = stpc;
        } else if (stp > stx)
            stpf = stpmax;
        else {
            stpf = stpmin;
        }
    }

    if (fp > fx) {
        sty = stp;
        fy = fp;
        dy = dp;
    } else {
        if (sgnd < 0.0) {
            sty = stx;
            fy = fx;
            dy = dx;
        }

        stx = stp;
        fx = fp;
        dx = dp;
    }

    stpf = std::min(stpmax, stpf);
    stpf = std::max(stpmin, stpf);
    stp = stpf;

    if (brackt & bound) {
        if (sty > stx) {
            stp = std::min(stx + 0.66 * (sty - stx), stp);
        } else {
            stp = std::max(stx + 0.66 * (sty - stx), stp);
        }
    }

    return 0;
}
} // namespace magnumafcpp
