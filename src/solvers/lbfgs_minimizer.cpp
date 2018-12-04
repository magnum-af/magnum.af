#include "lbfgs_minimizer.hpp"

LBFGS_Minimizer::LBFGS_Minimizer()
{
}

//Energy calculation
double LBFGS_Minimizer::E(const State& state){
  double solution = 0.;
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

af::array LBFGS_Minimizer::Gradient(const State& state){
    return cross4(state.m, cross4(state.m, Heff(state)));
    //TODO// this works!!!//return cross4(state.m, cross4(state.m, Heff(state)));
    //TODO//this is wrong, too slow maybe?//return - state.param.alpha*state.param.gamma/(1.+pow(state.param.alpha,2)) * cross4(state.m, cross4(state.m, Heff(state)));
}

double mydot (const af::array& a, const af::array& b){
    return full_inner_product(a, b);
}
double mynorm(const af::array &a) {
  return sqrt( mydot(a,a) );
}

double LBFGS_Minimizer::mxmxhMax(const State& state) {
    return maxnorm(Gradient(state));
}

/// LBFGS minimizer from Thomas Schrefl's bvec code
double LBFGS_Minimizer::Minimize(State& state){
    af::timer timer = af::timer::start();
    //af::print("h in minimize", af::mean(Heff(state)));//TODEL
    //af::print("Gradient", Gradient(state));//TODEL

    double eps  = 2.22e-16;
    double eps2 = sqrt(eps);
    double epsr = pow(eps,0.9);
    double tolf = 1e-6; //TODO find value of this->settings_.gradTol;
    double tolf2 = sqrt(tolf);
    double tolf3 = pow(tolf,0.3333333333333333333333333);
    double f = this->E(state);
    //double f = objFunc.both(x0, grad);// objFunc.both calcs Heff and E for not calculating Heff double
    // NOTE: objFunc.both calcs gradient and energy E
    af::array grad = this->Gradient(state);
    double gradNorm = maxnorm(grad);
    //af::print("grad", grad);
    std::cout << "f= " << f << std::endl;
    if ( gradNorm < (epsr*(1+fabs(f))) ) {
      return f;
    }
    const size_t m = 5;

    // array container
    af::array af_zero = (af::constant(0.0, state.mesh.dims, f64));
    std::array<af::array, m> sVector {{af_zero, af_zero, af_zero, af_zero, af_zero}};
    std::array<af::array, m> yVector {{af_zero, af_zero, af_zero, af_zero, af_zero}};

    // vector container
    //std::vector<af::array> sVector; //= af::constant(0.0, state.mesh.dims, f64);
    //std::vector<af::array> yVector; //= af::constant(0.0, state.mesh.dims, f64);
    //for (size_t i = 0; i < m; i++){
    //    sVector.push_back(af::constant(0.0, state.mesh.dims, f64));
    //    yVector.push_back(af::constant(0.0, state.mesh.dims, f64));
    //}

    std::array<double, m> alpha = {0.,0.,0.,0.,0.};
    af::array q = af::constant(0.0, state.mesh.dims, f64);
    af::array s = af::constant(0.0, state.mesh.dims, f64);
    af::array y = af::constant(0.0, state.mesh.dims, f64);
    double f_old = 0.0;
    af::array grad_old = af::constant(0.0, state.mesh.dims, f64);
    af::array x0 = state.m;
    af::array x_old = x0;

    size_t iter = 0, globIter = 0; 
    double H0k = 1;
        do {
            //TODO TODEL
            //
            //for (size_t i = 0; i < m; i++){
            //    printf("size svec = %i \n", sVector.size());
            //    af::print("svec", af::mean(sVector[i],0));
            //}
            //END TODO
            int cgSteps = 0;
            //this->minIterCount_++;
            f_old = f; 
            x_old = x0;
            //af::print("x0", af::mean(x0,0));
            grad_old = grad;

            q = grad;
            const int k = std::min(m, iter);
            for (int i = k - 1; i >= 0; i--) {
                const double rho = 1.0 / mydot(sVector[i],yVector[i]);
                alpha[i] = rho * mydot(sVector[i],q);
                q -= alpha[i] * yVector[i];
            }
            q = H0k * q; // NOTE: cg step skipped and only used else
            for (int i = 0; i < k; i++) {
                const double rho = 1.0 / mydot(sVector[i],yVector[i]);
                const double beta = rho * mydot(yVector[i],q);
                q += sVector[i] * (alpha[i] - beta);
            }

            // line 291 in .fe: decent = -grad.dot(q) 
            double phiPrime0 = -mydot(grad,q);
            if (phiPrime0 > 0) {
                q = grad;
                iter = 0;
                if (this->verbose_>2) {
                  std::cout << "descent " << std::endl;
                }
                phiPrime0 = -mydot(grad,q);
            }

            //TODO objFunc.updateTol(gradNorm);
            //TODO check version Schrefl vs Flo
            if ( -mydot(grad, q) > -1e-15){
               gradNorm = maxnorm(grad); 
               if (gradNorm < eps*(1.+fabs(f))){
                   std::cout << "Minimizer: Convergence reached (due to almost zero gradient (|g|=" << gradNorm << " < " << eps*(1.+fabs(f))<< ")!" << std::endl;
               }
            }
            //const double rate = MyMoreThuente<T, decltype(objFunc), 1>::linesearch(f, x_old, x0, grad, -q, objFunc, tolf);
            const double rate = linesearch(state, f, x_old, x0, grad, -q, tolf);
            //TODO/old/todel//const double rate = this->cvsrch(state, x_old, x0, f, grad, -q, tolf);
            std::cout.precision(24);
            std::cout << "rate= " << rate << std::endl;
            state.m = x0;
            if (this->verbose_>3 && rate == 0.0) {
              std::cout << bold_red("Error: LBFGS_Minimizer: linesearch failed") << std::endl;
              exit(0);
            }
            double f1 = 1+fabs(f);
            double gradNorm = maxnorm(grad);
            if (gradNorm < (epsr*f1)) {
              return f;
            }
            s = x0 - x_old;
            if (this->verbose_ > 1) {
              //std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
              std::cout << "bb> " << globIter << " " << f << " " << " " << gradNorm << " " << cgSteps << " " << rate << std::endl;
            }
            if ( ( (f_old-f)   < (tolf*f1) ) && 
                 (  maxnorm(s) < (tolf2*(1+maxnorm(x0))) ) && 
                 (  gradNorm  <= (tolf3*f1) ) )  {
              break;
            }
            y = grad - grad_old;
            double ys = mydot(y,s);
            //TODO this is always true, ( because rate = 0? )
            if (ys <= eps2*mynorm(y)*mynorm(s)) { // Dennis, Schnabel 9.4.1
              if (this->verbose_>2) {
                std::cout << iter << " skip update!!!!!!!!! " << std::endl;
              }
            }
            else { 
              if (iter < m) {
                //af::print("s", af::mean(s,0));
                sVector[iter] = s;
                //af::print("svec", af::mean(sVector[iter],0));
                yVector[iter] = y;
              } else {
                std::rotate(sVector.begin(),sVector.begin()+1,sVector.end());
                sVector[m-1] = s;
                std::rotate(yVector.begin(),yVector.begin()+1,yVector.end());
                yVector[m-1] = y;
              }
              H0k = ys / mydot(y,y);
              iter++;
            }

            globIter++;

        } while (globIter < this->maxIter_); 
        if (globIter >= this->maxIter_ && this->maxIter_ > 99) {
          std::cout << "WARNING : maximum number of iterations exceeded in LBFGS" << std::endl;
        }
        if (this->verbose_ > 0) {
            State temp = state;
            temp.m = grad;
          auto mxh = this->mxmxhMax(temp);
          auto deltaF = f-f_old;
          std::cout << "aa> " << globIter << " " << deltaF << " " << " " << mxh << std::endl;
        }
      return f;        



    std::cout << "LBFGS_Minimizer: minimize in [s]: " << af::timer::stop(timer) << std::endl;
}; 

double LBFGS_Minimizer::linesearch(const State& state, double &fval, const af::array &x_old, af::array &x, af::array &g, const af::array &searchDir, double tolf) { 
    double ak = 1.0;
    //TODO//this in not changing, (because rate becomes zero an nothing changes//af::print("x", af::mean(x,0));
    //TODEL//af::print("x_old", af::mean(x_old,0)); 
    //TODEL//af::print("x", af::mean(x,0)); 
    //TODEL//af::print("g", af::mean(g,0)); 
    //TODEL//af::print("searchDir", af::mean(searchDir,0)); 
    //TODEL//std::cout << "fval= " << fval << " ak= " << ak << std::endl;
    cvsrch(state, x_old, x, fval, g, ak, searchDir, tolf);
    //TODEL//af::print("x_old", af::mean(x_old,0)); 
    //TODEL//af::print("x", af::mean(x,0)); 
    //TODEL//af::print("g", af::mean(g,0)); 
    //TODEL//af::print("searchDir", af::mean(searchDir,0)); 
    //TODEL//std::cout << "fval= " << fval << " ak= " << ak << std::endl;
    //TODO//this in not changing, (because rate becomes zero an nothing changes//af::print("x", af::mean(x,0));
    return ak;
}

//static int cvsrch(P &objFunc, const vex::vector<double> &wa, vex::vector<double> &x, double &f, vex::vector<double> &g, double &stp, const vex::vector<double> &s, double tolf) {
//TODO//TODEL//int LBFGS_Minimizer::cvsrch(const State& state, const af::array &wa, af::array &x, double &f, af::array &g, const af::array &s, double tolf) {// ak = 1.0 == double &stp moved into function
int LBFGS_Minimizer::cvsrch(const State& state, const af::array &wa, af::array &x, double &f, af::array &g, double &stp, const af::array &s, double tolf) {
  // we rewrite this from MIN-LAPACK and some MATLAB code

  int info           = 0;
  int infoc          = 1;

//  const int n        = x.dims(0)*x.dims(1)*x.dims(2)*x.dims(3); (void) n;
  const double xtol   = 1e-15;
  const double ftol   = 1.0e-4;  // c1
  const double gtol   = 0.9;   // c2
  const double eps    = tolf;
  const double stpmin = 1e-15;
  const double stpmax = 1e15;
  const double xtrapf = 4;
  const int maxfev   = 20;
  int nfev           = 0;

  double dginit = mydot(g,s);
  if (dginit >= 0.0) {
    // no descent direction
    // TODO: handle this case
    std::cout << " no descent " << dginit << std::endl;
    return -1;
  }

  bool brackt      = false;
  bool stage1      = true;

  double finit      = f;
  double dgtest     = ftol * dginit;
  double width      = stpmax - stpmin;
  double width1     = 2 * width;
  // vex::vector<double> wa(x.queue_list(),x.size());
  // wa = x;

  double stx        = 0.0;
  double fx         = finit;
  double dgx        = dginit;
  double sty        = 0.0;
  double fy         = finit;
  double dgy        = dginit;

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
    if (
    (brackt && ((stp <= stmin) | (stp >= stmax)))
    | (nfev >= maxfev - 1 ) | (infoc == 0)
    | (brackt & (stmax - stmin <= xtol * stmax))) {
      stp = stx;
    }

    //// test new point
    //// x = wa + stp * s;
    //// objFunc.normalizeVector(x);
    //objFunc.update(stp,wa,s,x);
    //TODO check
    x = wa + stp * s; // this is equivalent to objFunc.update(stp,wa,s,x);
    //TODO/TODEL//af::print("linesearch x", af::mean(x,0));
    x = renormalize_handle_zero_values(x);
    //TODO/TODEL//af::print("linesearch x renorm", af::mean(x,0));

    // f = objFunc.both(x, g);
    // NOTE: objFunc.both calcs gradient and energy E
    State state_x = state;
    state_x.m = x;
    f = this->E(state_x);
    //TODEL//This is changing//std::cout<< "in linesearch f= " << f << std::endl;
    g = this->Gradient(state_x);
    //TODO//TODEL//this ich changing//af::print("g", af::mean(g,0));
    //END TODO check
    nfev++;
    double dg = mydot(g,s);
    double ftest1 = finit + stp * dgtest;
    double ftest2 = finit + eps*fabs(finit);
    double ft = 2*ftol-1;

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
      double fm = f - stp * dgtest;
      double fxm = fx - stx * dgtest;
      double fym = fy - sty * dgtest;
      double dgm = dg - dgtest;
      double dgxm = dgx - dgtest;
      double dgym = dgy - dgtest;

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

  std::cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXX why I am here ? XXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxx " << std::endl;
  exit(0);
  return 0;
}

int LBFGS_Minimizer::cstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, double& stp,
double& fp, double& dp, bool& brackt, double& stpmin, double& stpmax, int& info) {
  info = 0;
  bool bound = false;

  // Check the input parameters for errors.
  if ((brackt & ((stp <= std::min(stx, sty) ) | (stp >= std::max(stx, sty)))) | (dx * (stp - stx) >= 0.0)
  | (stpmax < stpmin)) {
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
    double gamma = s * sqrt((theta / s) * (theta / s)  - (dx / s) * (dp / s));
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
    double s = std::max(theta, std::max( dx, dp));
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
