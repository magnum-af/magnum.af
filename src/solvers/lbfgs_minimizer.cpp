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
        std::cout<<"ERROR: LBFGS_Minimizer::Heff: Number of _llgterms == 0. Please add at least one term to LBFGS_Minimizer.llgterms_! Aborting..."<<std::endl;
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
    return - state.param.alpha*state.param.gamma/(1.+pow(state.param.alpha,2)) * cross4(state.m, cross4(state.m, Heff(state)));
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
    double tolf = 1e-12; //TODO find value of this->settings_.gradTol;
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

    std::vector<af::array> sVector; //= af::constant(0.0, state.mesh.dims, f64);
    std::vector<af::array> yVector; //= af::constant(0.0, state.mesh.dims, f64);
    for (size_t i = 0; i < m; i++){
        sVector.push_back(af::constant(0.0, state.mesh.dims, f64));
        yVector.push_back(af::constant(0.0, state.mesh.dims, f64));
    }

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
            int cgSteps = 0;
            //this->minIterCount_++;
            f_old = f; 
            x_old = x0;
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

            // Decent = -grad.dot(q) in .fe
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
            //const double rate = MyMoreThuente<T, decltype(objFunc), 1>::linesearch(f, x_old, x0, grad, -q, objFunc, tolf);
            const double rate = 1.0; // TODO temp-fix proposed by Flo 
            if (this->verbose_>3 && rate == 0.0) {
              std::cout << "linesearch failed" << std::endl;
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
            if (ys <= eps2*mynorm(y)*mynorm(s)) { // Dennis, Schnabel 9.4.1
              if (this->verbose_>2) {
                std::cout << iter << " skip update!!!!!!!!! " << std::endl;
              }
            }
            else { 
              if (iter < m) {
                sVector[iter] = s;
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
