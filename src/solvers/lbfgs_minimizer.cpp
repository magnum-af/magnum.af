#include "lbfgs_minimizer.hpp"
#include "../state.hpp"
#include "../misc.hpp"
#include "../func.hpp"
#include "../llg_terms/LLGTerm.hpp"
#include "arrayfire.h"
#include <iomanip>
#include <memory>
#include <list>
#include <algorithm>

namespace magnumaf{


LBFGS_Minimizer::LBFGS_Minimizer(float tolerance, size_t maxIter, int verbose): tolerance_(tolerance), maxIter_(maxIter), verbose(verbose)
{
}

LBFGS_Minimizer::LBFGS_Minimizer(LlgTerms llgterms, float tolerance, size_t maxIter, int verbose): llgterms_(llgterms), tolerance_(tolerance), maxIter_(maxIter), verbose(verbose)
{
}


//Energy calculation
float LBFGS_Minimizer::E(const State& state){
  float solution = 0.;
  for(unsigned i=0;i<llgterms_.size();++i){
    solution+=llgterms_[i]->E(state);
  }
  return solution;
}
// Calculation of effective field
af::array LBFGS_Minimizer::Heff(const State& state){
    if( llgterms_.size() == 0){
        std::cout<< bold_red("ERROR: LBFGS_Minimizer::Heff: Number of _llgterms == 0. Please add at least one term to LBFGS_Minimizer.llgterms_! Aborting...")<<std::endl;
        exit (EXIT_FAILURE);
     }
    af::timer timer=af::timer::start();
    af::array solution = llgterms_[0]->h(state);
    for(unsigned i = 1; i < llgterms_.size() ; ++i ){
        solution+=llgterms_[i]->h(state);
    }
    time_calc_heff_ += af::timer::stop(timer);
    return solution;
}

float LBFGS_Minimizer::EnergyAndGradient(const State& state, af::array& gradient){
    if( llgterms_.size() == 0){
        std::cout<< bold_red("ERROR: LBFGS_Minimizer::Heff: Number of _llgterms == 0. Please add at least one term to LBFGS_Minimizer.llgterms_! Aborting...")<<std::endl;
        exit (EXIT_FAILURE);
     }
    af::timer timer=af::timer::start();
    //Avoiding array with zeros, starting loop with second term in llgterms
    af::array h = llgterms_[0]->h(state);
    float energy = llgterms_[0]->E(state, h);
    for(unsigned i = 1; i < llgterms_.size() ; ++i ){
    //Alternative:
    //af::array h = af::constant(0, state.mesh.dims, f32);//llgterms_[0]->h(state);
    //float energy = 0;//llgterms_[0]->E(state, h);
    //for(unsigned i = 0; i < llgterms_.size() ; ++i ){
        af::array temp_h = llgterms_[i]->h(state);
        h+=temp_h;
        energy+=llgterms_[i]->E(state, temp_h);
    }
    gradient = 1./(constants::mu0 * state.Ms) * cross4(state.m, cross4(state.m, h));
    time_calc_heff_ += af::timer::stop(timer);
    return energy;
}

af::array LBFGS_Minimizer::Gradient(const State& state){
    //TODO//this runs, check correct way//return cross4(state.m, cross4(state.m, Heff(state)));
    return 1./(constants::mu0 * state.Ms) * cross4(state.m, cross4(state.m, Heff(state)));// TODO elaborate correct way

    //return cross4(state.m, cross4(state.m, constants::mu0 * state.Ms * state.mesh.V * Heff(state)));
    //return constants::mu0 * state.Ms * state.mesh.V * cross4(state.m, cross4(state.m, Heff(state)));
    //TODO//this runs//return cross4(state.m, cross4(state.m, Heff(state)));
    //TODO// this works!!!//return cross4(state.m, cross4(state.m, Heff(state)));
    //TODO//this is wrong, too slow maybe?//
    //return pow(constants::mu0, 2)/state.mesh.V * state.Ms * cross4(state.m, cross4(state.m, Heff(state)));
    //return state.Ms * cross4(state.m, cross4(state.m, Heff(state)));
}

float mydot (const af::array& a, const af::array& b){
    return full_inner_product(a, b);
}
float mynorm(const af::array &a) {
  return sqrt( mydot(a, a) );
}

float LBFGS_Minimizer::mxmxhMax(const State& state) {
    return maxnorm(Gradient(state));
}

/// LBFGS minimizer from Thomas Schrefl's bvec code
//TODO Currently Minimize() fails when called second time, i.e. when m already relaxed linesearch rate is 0 even with changed zeeman field
float LBFGS_Minimizer::Minimize(State& state){
    std::cout.precision(24);
    af::timer timer = af::timer::start();
    //af::print("h in minimize", af::mean(Heff(state)));//TODEL
    //af::print("Gradient", Gradient(state));//TODEL

    float eps  = 2.22e-16;
    float eps2 = sqrt(eps);
    float epsr = pow(eps, 0.9);
    //float tolerance_ = 1e-6; //TODO find value of this->settings_.gradTol;
    float tolf2 = sqrt(tolerance_);
    float tolf3 = pow(tolerance_, 0.3333333333333333333333333);
    //float f = this->E(state);
    //af::array grad = this->Gradient(state);
    af::array grad;
    float f = this->EnergyAndGradient(state, grad);
    //float f = objFunc.both(x0, grad);// objFunc.both calcs Heff and E for not calculating Heff float
    // NOTE: objFunc.both calcs gradient and energy E
    float gradNorm = maxnorm(grad);
    //af::print("grad", grad);
    if(this->verbose > 0){std::cout << "f= " << f << std::endl;}
    if ( gradNorm < (epsr*(1+fabs(f))) ) {
      return f;
    }
    const size_t m = 5;

    // array container
    af::array af_zero = (af::constant(0.0, state.mesh.dims, f32));
    std::array<af::array, m> sVector {{af_zero, af_zero, af_zero, af_zero, af_zero}};
    std::array<af::array, m> yVector {{af_zero, af_zero, af_zero, af_zero, af_zero}};

    // vector container
    //std::vector<af::array> sVector; //= af::constant(0.0, state.mesh.dims, f32);
    //std::vector<af::array> yVector; //= af::constant(0.0, state.mesh.dims, f32);
    //for (size_t i = 0; i < m; i++){
    //    sVector.push_back(af::constant(0.0, state.mesh.dims, f32));
    //    yVector.push_back(af::constant(0.0, state.mesh.dims, f32));
    //}

    std::array<float, m> alpha = {0., 0., 0., 0., 0.};
    af::array q = af::constant(0.0, state.mesh.dims, f32);
    af::array s = af::constant(0.0, state.mesh.dims, f32);
    af::array y = af::constant(0.0, state.mesh.dims, f32);
    float f_old = 0.0;
    af::array grad_old = af::constant(0.0, state.mesh.dims, f32);
    af::array x_old = state.m;

    size_t iter = 0, globIter = 0;
    float H0k = 1;
        do {
            //TODO TODEL
            //
            //for (size_t i = 0; i < m; i++){
            //    printf("size svec = %i \n", sVector.size());
            //    af::print("svec", af::mean(sVector[i], 0));
            //}
            //END TODO
            int cgSteps = 0;
            //this->minIterCount_++;
            f_old = f;
            x_old = state.m;
            //af::print("state.m", af::mean(state.m, 0));
            grad_old = grad;

            q = grad;
            const int k = std::min(m, iter);
            for (int i = k - 1; i >= 0; i--) {
                const float rho = 1.0 / mydot(sVector[i], yVector[i]);
                alpha[i] = rho * mydot(sVector[i], q);
                q -= alpha[i] * yVector[i];
            }
            q = H0k * q; // NOTE: cg step skipped and only used else
            for (int i = 0; i < k; i++) {
                const float rho = 1.0 / mydot(sVector[i], yVector[i]);
                const float beta = rho * mydot(yVector[i], q);
                q += sVector[i] * (alpha[i] - beta);
            }

            // line 291 in .fe: decent = -grad.dot(q)
            float phiPrime0 = -mydot(grad, q);
            if (phiPrime0 > 0) {
                q = grad;
                iter = 0;
                if (this->verbose>2) {
                  std::cout << "descent " << std::endl;
                }
                phiPrime0 = -mydot(grad, q);
            }

            //TODO objFunc.updateTol(gradNorm);
            //TODO check version Schrefl vs Flo
            if ( -mydot(grad, q) > -1e-15){
               gradNorm = maxnorm(grad);
               if (gradNorm < eps*(1.+fabs(f)) and this->verbose > 0){
                   std::cout << "Minimizer: Convergence reached (due to almost zero gradient (|g|=" << gradNorm << " < " << eps*(1.+fabs(f))<< ")!" << std::endl;
               }
            }
            //const float rate = MyMoreThuente<T, decltype(objFunc), 1>::linesearch(f, x_old, x0, grad, -q, objFunc, tolf);
            const float rate = linesearch(state, f, x_old, grad, -q, tolerance_);
            //TODO/old/todel//const float rate = this->cvsrch(state, x_old, x0, f, grad, -q, tolf);
            if (rate == 0.0) {
                if(this->verbose>0 ){
                    std::cout << red("Warning: LBFGS_Minimizer: linesearch returned rate == 0.0. This should not happen, elaborate! (maybe m is already relaxed?)") << std::endl;
                    //std::cout << bold_red("Error: LBFGS_Minimizer: linesearch failed, rate == 0.0") << std::endl;
                }
              //TODO//exit(0);// why exit with 0 (= sucess)? maybe change to return(0) or, (if necessary)//exit(EXIT_FAILURE);
            }
            float f1 = 1+fabs(f);
            float gradNorm = maxnorm(grad);
            if (gradNorm < (epsr*f1)) {
              return f;
            }
            s = state.m - x_old;
            if (this->verbose > 1) {
              //std::cout << std::setprecision(std::numeric_limits<float>::digits10 + 1);
              std::cout << std::setprecision(std::numeric_limits<float>::digits10+12);
              std::cout << "bb> " << globIter << " " << f << " " << " " << gradNorm << " " << cgSteps << " " << rate << std::endl;
            }
            if (of_convergence.is_open()){
                of_convergence << (f_old-f)/(tolerance_*f1)<< "\t" << maxnorm(s) / (tolf2*(1+maxnorm(state.m)))<< "\t" << gradNorm/(tolf3*f1) << std::endl;
            }
            if ( ( (f_old-f)   < (tolerance_*f1) ) &&
                 (  maxnorm(s) < (tolf2*(1+maxnorm(state.m))) ) &&
                 (  gradNorm  <= (tolf3*f1) ) )  {
              break;
            }
            y = grad - grad_old;
            float ys = mydot(y, s);
            if (ys <= eps2*mynorm(y)*mynorm(s)) { // Dennis, Schnabel 9.4.1
              if (this->verbose>2) {
                std::cout << iter << red("WARNING: LBFGS_Minimizer:: skipping update!") << std::endl;
              }
            }
            else {
              if (iter < m) {
                //af::print("s", af::mean(s, 0));
                sVector[iter] = s;
                //af::print("svec", af::mean(sVector[iter], 0));
                yVector[iter] = y;
              } else {
                std::rotate(sVector.begin(), sVector.begin()+1, sVector.end());
                sVector[m-1] = s;
                std::rotate(yVector.begin(), yVector.begin()+1, yVector.end());
                yVector[m-1] = y;
              }
              H0k = ys / mydot(y, y);
              iter++;
            }

            globIter++;

        } while (globIter < this->maxIter_);
        if (globIter >= this->maxIter_ && this->maxIter_ > 99 and this->verbose > 0) {
          std::cout << "WARNING : maximum number of iterations exceeded in LBFGS" << std::endl;
        }
        if (this->verbose > 0) {
            State temp = state;
            temp.m = grad;
          auto mxh = this->mxmxhMax(temp);
          auto deltaF = f-f_old;
          std::cout << "aa> " << globIter << " " << deltaF << " " << " " << mxh << std::endl;
        }
      return f;



    if( this->verbose > 0) {std::cout << "LBFGS_Minimizer: minimize in [s]: " << af::timer::stop(timer) << std::endl;}
};

float LBFGS_Minimizer::linesearch(State& state, float &fval, const af::array &x_old, af::array &g, const af::array &searchDir, const float tolf) {
    float rate = 1.0;
    cvsrch(state, x_old, fval, g, rate, searchDir, tolf);
    return rate;
}

//static int cvsrch(P &objFunc, const vex::vector<float> &wa, vex::vector<float> &x, float &f, vex::vector<float> &g, float &stp, const vex::vector<float> &s, float tolf) {
//TODO//TODEL//int LBFGS_Minimizer::cvsrch(const State& state, const af::array &wa, af::array &x, float &f, af::array &g, const af::array &s, float tolf) {// ak = 1.0 == float &stp moved into function
int LBFGS_Minimizer::cvsrch(State& state, const af::array &wa, float &f, af::array &g, float &stp, const af::array &s, const float tolf) {
  // we rewrite this from MIN-LAPACK and some MATLAB code

  int info           = 0;
  int infoc          = 1;

//  const int n        = x.dims(0)*x.dims(1)*x.dims(2)*x.dims(3); (void) n;
  const float xtol   = 1e-15;
  const float ftol   = 1.0e-4;  // c1
  const float gtol   = 0.9;   // c2
  const float eps    = tolf;
  const float stpmin = 1e-15;
  const float stpmax = 1e15;
  const float xtrapf = 4;
  const int maxfev   = 20;
  int nfev           = 0;

  float dginit = mydot(g, s);
  if (dginit >= 0.0) {
    // no descent direction
    // TODO: handle this case
    std::cout << red("WARNING: LBFGS_Minimizer:: no descent ") << dginit << std::endl;
    return -1;
  }

  bool brackt      = false;
  bool stage1      = true;

  float finit      = f;
  float dgtest     = ftol * dginit;
  float width      = stpmax - stpmin;
  float width1     = 2 * width;
  // vex::vector<float> wa(x.queue_list(), x.size());
  // wa = x;

  float stx        = 0.0;
  float fx         = finit;
  float dgx        = dginit;
  float sty        = 0.0;
  float fy         = finit;
  float dgy        = dginit;

  float stmin;
  float stmax;

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
    if (
    (brackt && ((stp <= stmin) | (stp >= stmax)))
    | (nfev >= maxfev - 1 ) | (infoc == 0)
    | (brackt & (stmax - stmin <= xtol * stmax))) {
        stp = stx;
        if( this->verbose > 0) std::cout << "NOTE: LBFGS_Minimizer:: Oops, stp= " << stp << std::endl;
    }

    //// test new point
    //// x = wa + stp * s;
    //// objFunc.normalizeVector(x);
    //objFunc.update(stp, wa, s, x);
    state.m = wa + stp * s; //TODO check// this should be equivalent to objFunc.update(stp, wa, s, x);
    state.m = renormalize_handle_zero_values(state.m);
    f = this->EnergyAndGradient(state, g);
    nfev++;
    float dg = mydot(g, s);
    float ftest1 = finit + stp * dgtest;
    float ftest2 = finit + eps*fabs(finit);
    float ft = 2*ftol-1;

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

    if (((f<=ftest2)&&(ft*dginit>=dg)) && (fabs(dg) <= gtol * (-dginit)))
      info = 1;

    // terminate when convergence reached
    if (info != 0) {
      return -1;
    }

    // if (stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol)*dginit))
    if (stage1 & ((f<=ftest2)&&(ft*dginit>=dg)) & (dg >= std::min(ftol, gtol)*dginit))  // approx wolfe 1
      stage1 = false;

    // if (stage1 & (f <= fx) & (f > ftest1)) {
    if (stage1 & (f <= fx) & (not ((f<=ftest2)&&(ft*dginit>=dg)))) { // not wolfe 1 --> not approx wolfe 1
      float fm = f - stp * dgtest;
      float fxm = fx - stx * dgtest;
      float fym = fy - sty * dgtest;
      float dgm = dg - dgtest;
      float dgxm = dgx - dgtest;
      float dgym = dgy - dgtest;

      int rsl = cstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc); (void) rsl;

      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    } else {
      // this is ugly and some variables should be moved to the class scope
      int rsl = cstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc); (void) rsl;
    }

    if (brackt) {
      if (fabs(sty - stx) >= 0.66 * width1)
        stp = stx + 0.5 * (sty - stx);
      width1 = width;
      width = fabs(sty - stx);
    }
  }

  std::cout << bold_red("ERROR: LBFGS_Minimizer: cvsrch: XXXXXXXXXXXXXXXXXXXXXXXXXXXX why I am here ? XXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxx ") << std::endl;
  exit(0);
  return 0;
}

int LBFGS_Minimizer::cstep(float& stx, float& fx, float& dx, float& sty, float& fy, float& dy, float& stp,
float& fp, float& dp, bool& brackt, float& stpmin, float& stpmax, int& info) {
  info = 0;
  bool bound = false;

  // Check the input parameters for errors.
  if ((brackt & ((stp <= std::min(stx, sty) ) | (stp >= std::max(stx, sty)))) | (dx * (stp - stx) >= 0.0)
  | (stpmax < stpmin)) {
    return -1;
  }

  float sgnd = dp * (dx / fabs(dx));

  float stpf = 0;
  float stpc = 0;
  float stpq = 0;

  if (fp > fx) {
    info = 1;
    bound = true;
    float theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
    float s = std::max(theta, std::max(dx, dp));
    float gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
    if (stp < stx)
      gamma = -gamma;
    float p = (gamma - dx) + theta;
    float q = ((gamma - dx) + gamma) + dp;
    float r = p / q;
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
    float theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    float s = std::max(theta, std::max(dx, dp));
    float gamma = s * sqrt((theta / s) * (theta / s)  - (dx / s) * (dp / s));
    if (stp > stx)
      gamma = -gamma;

    float p = (gamma - dp) + theta;
    float q = ((gamma - dp) + gamma) + dx;
    float r = p / q;
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
    float theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
    float s = std::max(theta, std::max( dx, dp));
    float gamma = s * sqrt(std::max((float) 0., (theta / s) * (theta / s) - (dx / s) * (dp / s)));
    if (stp > stx)
      gamma = -gamma;
    float p = (gamma - dp) + theta;
    float q = (gamma + (dx - dp)) + gamma;
    float r = p / q;
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
      float theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
      float s = std::max(theta, std::max(dy, dp));
      float gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
      if (stp > sty)
        gamma = -gamma;

      float p = (gamma - dp) + theta;
      float q = ((gamma - dp) + gamma) + dy;
      float r = p / q;
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
      stp = std::min(stx + (float) 0.66 * (sty - stx), stp);
    } else {
      stp = std::max(stx + (float) 0.66 * (sty - stx), stp);
    }
  }

  return 0;

}
}// namespace magnumaf
