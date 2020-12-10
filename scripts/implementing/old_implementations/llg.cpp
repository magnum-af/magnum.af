#include "llg.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>

using namespace magnumafcpp;

using namespace af;

// Energy calculation
double LLG::E(const State& state) {
    double solution = 0.;
    for (unsigned i = 0; i < Fieldterms.size(); ++i) {
        solution += Fieldterms[i]->Energy_in_J(state);
    }
    return solution;
}

// array LLG::givescale(const array& a){
//  return atol+rtol*abs(a);
//}
void LLG::write_fieldterms_atom(const State& state, const std::string filepath) {
    for (unsigned i = 0; i < Fieldterms.size(); ++i) {
        vti_writer_atom(Fieldterms[i]->H_in_Apm(state), state.mesh, filepath + std::to_string(i));
    }
}
void LLG::write_fieldterms_micro(const State& state, const std::string filepath) {
    for (unsigned i = 0; i < Fieldterms.size(); ++i) {
        vti_writer_micro(Fieldterms[i]->H_in_Apm(state), state.mesh, filepath + std::to_string(i));
    }
}

void LLG::relax(State& state, const double precision, const int iloop, const int iwritecout) {
    timer t = af::timer::start();
    double E_prev = 1e20;
    while (fabs((E_prev - E(state)) / E_prev) > precision) {
        E_prev = E(state);
        for (int i = 0; i < iloop; i++) {
            state.m = step(state);
        }
        if (state.steps % iwritecout == 0)
            std::cout << "step " << state.steps << " rdiff= " << fabs((E_prev - E(state)) / E_prev) << std::endl;
    }
    std::cout << "timerelax [af-s]: " << af::timer::stop(t) << ", current llg steps = " << state.steps << std::endl;
}

array LLG::fdmdt(const array& m, const array& heff) {
    calls_fdmdt++;
    timer_fdmdt = timer::start();
    if (fdmdt_dissipation_term_only) {
        dmdt = -state0.material.alpha * state0.material.gamma / (1. + pow(state0.material.alpha, 2)) *
               cross4(m, cross4(m, heff));
    } else {
        dmdt = -state0.material.gamma / (1. + pow(state0.material.alpha, 2)) * cross4(m, heff) -
               state0.material.alpha * state0.material.gamma / (1. + pow(state0.material.alpha, 2)) *
                   cross4(m, cross4(m, heff));
    }
    time_fdmdt += af::timer::stop(timer_fdmdt);
    return dmdt;
}

array LLG::fdmdtminimal(array m, array heff) { // array LLG::fdmdt(array& m, array& heff){
    dmdt = cross4(m, heff);
    return -state0.material.gamma / (1. + pow(state0.material.alpha, 2)) * cross4(m, heff) -
           state0.material.alpha * state0.material.gamma / (1. + pow(state0.material.alpha, 2)) * cross4(m, dmdt);
}

// Calculation of effective field
array LLG::fheff(const array& m) {
    array solution = constant(0., state0.mesh::dims_v(mesh), f64);
    timer_heff = timer::start();

    state0.m = m; // TODO avoid state0 in the first place
    // TODO avoid this line  State temp(state0.mesh, state0.material, m);
    for (unsigned i = 0; i < Fieldterms.size(); ++i) {
        solution += Fieldterms[i]->H_in_Apm(state0);
    }

    time_heff += af::timer::stop(timer_heff);
    return solution;
}

double LLG::cpu_time() {
    double cpu_time = 0.;
    for (unsigned i = 0; i < Fieldterms.size(); ++i) {
        cpu_time += Fieldterms[i]->elapsed_eval_time();
    }
    return cpu_time;
}
void LLG::print_cpu_time(std::ostream& stream) { stream << "cpu_time = " << LLG::cpu_time() << " [s]" << std::endl; }

long int LLG::h_addr(const State& state) {
    // array fheff_pass_to_python = fheff(state.m);
    // fheff_pass_to_python.lock(); //Caution with locking arrays
    // std::cout << "LLG::h_addr:"<< (long int) fheff_pass_to_python.get() <<
    // std::endl; return (long int) fheff_pass_to_python.get();

    //  h_addr_temp_array = fheff(state.m);
    //  //h_addr_temp_array.lock(); //check if it works without locking
    //  return (long int) h_addr_temp_array.get();

    // With vector, but also with sagfaults for second call
    h_addr_temp_array.push_back(fheff(state.m));
    h_addr_temp_array.back().lock();
    return (long int)h_addr_temp_array.back().get();
}

//  std::cout<<"cpu_time = "<<llg.cpu_time()<<""<<std::endl;
//// Calculation of effective field
// array LLG::fheffminimal(array m){
//  if (state0.b_zee){return Demag.h(m) + Exch.h(m) + state0.h_zee;}
//  return Demag.h(m) + Exch.h(m);
//}

//// Calculation of effective field with zeeman field
// array LLG::fheff(array m, array h_zee){
//  return Demag.solve(m) + Exch.solve(m) + h_zee;
//}

// bool LLG::controller(const double err, double& h){
//  static const double beta =0.4/5.0;
//  static const double alpha = 0.2 - 0.75*beta;
//  static const double headroom{0.9};
//  static const double minscale{0.2};
//  static const double maxscale{10.};
//
//  double scale;
//
////  std::cout << "err = " << err << " h= " << h << "reject "<< reject << "
/// t="<<true << " f=" <<false << std::endl;
//  if (err <= 1.0){
//    if (err == 0.0)
//      scale = maxscale;
//    else{
//      scale=headroom*pow(err, -alpha)*pow(errold, beta);
//      if (scale < minscale){
//        scale=minscale;
//        counter_maxscale++;
//      }
//      if (scale > maxscale){
//        scale=maxscale;
//        counter_minscale++;
//      }
//    }
//  if (reject)
//    hnext=h*std::min(scale, 1.0);//Do not increase stepsize if previous try
//    was rejected
//  else
//    hnext=h*scale;
//  if(hnext<=hmin) {
//    hnext=hmin;
//    counter_hmin++;
//    std::cout << "Warning: hnext reached hmin in if, error bounds may be
//    invalid"<<std::endl;
//  }
//  if(hnext>=hmax){
//    hnext=hmax;
//    counter_hmax++;
//    std::cout << "Warning: hnext reached hmax in if, error bounds may be
//    invalid"<<std::endl;
//  }
//  //if(hnext<=hmin) hnext=hmin, std::cout << "Warning: hmin reached in if,
//  error bounds may be invalid"<<std::endl;
//  //if(hnext>=hmax) hnext=hmax, std::cout << "Warning: hmax reached in if,
//  error bounds may be invalid"<<std::endl;
//
//  errold=std::max(err, 1.0e-4);//Why?
//  reject=false;
//  counter_accepted++;
//  return true;
//  }
//  else{
//    scale=std::max(headroom*pow(err, -alpha), minscale);
//    h *= scale;
//    if(h<=hmin) {
//      h=hmin;
//      counter_hmin++;
//      std::cout << "Warning: hmin reached in else, error bounds may be
//      invalid"<<std::endl; return true;
//    }
//    if(h>=hmax){
//      h=hmax;
//      counter_hmax++;
//      std::cout << "Warning: hmax reached in else, error bounds may be
//      invalid"<<std::endl; return true;
//    }
//    reject=true;
//    counter_reject++;
//    return false;
//  }
//}

array LLG::step(State& state) {
    timer_integrator = timer::start();
    // h=1e-12;
    // hnext=h;
    // Fixed Stepsize Integrators
    if (state.material.mode == 0) {
        mtemp = explicitEuler(state.m, h);
    } else if (state.material.mode == 1) {
        mtemp = rk4(state.m, h);
    } else if (state.material.mode == 2) {
        mtemp = rk4minimal(state.m, h);
    } else if (state.material.mode == 3) {
        mtemp = rk4_3o8(state.m, h);
    } else if (state.material.mode == 5) {
        mtemp = RKF5(state.m, h);
    }
    // Adaptive Stepsize Integrators
    else {
        int while_break = 0;
        do {
            while_break++;
            if (while_break > 100)
                std::cout << "Warning: While_break > 100, break called" << std::endl;
            // BS3
            if (state.material.mode == 4) {
                mtemp = BS23(state.m, h, err);
            }
            // RKF
            else if (state.material.mode == 6) {
                mtemp = RKF45(state.m, h, err);
            }
            // CK5
            else if (state.material.mode == 7) {
                mtemp = CK45(state.m, h, err);
            }
            // Tsit45
            else if (state.material.mode == 8) {
                mtemp = tsit45(state.m, h, err);
            }
            // DP5
            else if (state.material.mode == 9) {
                // rk_rel_tol_error=5.e-2;
                mtemp = DP45(state.m, h, err);
            }
            // BS5
            else if (state.material.mode == 10) {
                mtemp = BS45(state.m, h, err);
            }
            // BS45 double error
            else if (state.material.mode == 11) {
                mtemp = BS45de(state.m, h, err);
            }
            // DP8
            else if (state.material.mode == 12) {
                mtemp = DP78(state.m, h, err);
            }
            if (controller.success(err, h))
                break;

        } while (while_break < 100);
    }
    state.t += h;
    h_stepped = h;
    state0.t += h; // TODO avoid state0
    mtemp += state.m;
    h = controller.get_hnext();
    if (state0.state.material.afsync)
        af::sync();
    time_integrator += timer::stop(timer_integrator);
    // mean(mtemp, 3); //TODO this is needed to avoid cuda crash?!?
    state.steps++;
    calls++;

    // Normalization
    // return normalize(mtemp);
    ////TODO better handle?
    if (state.Ms.isempty())
        return normalize(mtemp);
    else
        return (normalize_handle_zero_vectors(mtemp)); // Normalized array where all initial values == 0 are set
                                                       // to 0
    // return mtemp;

    // TODO
    //    if(fabs(maxnorm(vecnorm(mtemp))-1.) > llg_normtol){
    //      llg_normalize_counter++;
    //      llg_wasnormalized=true;
    //      std::cout<<"Normalization: "<<maxnorm(vecnorm(mtemp))<<"
    //      Counter="<<llg_normalize_counter<<std::endl; return
    //      normalize(mtemp);
    //    }
    //    else{
    //      llg_wasnormalized=false;
    //      return mtemp;
    //    }
};

// Time Integrators
array LLG::explicitEuler(const array& m, double dt) {
    heff = fheff(m);
    return dt * fdmdt(m, heff);
}

// Bogacki-Shampine 2/3rd order  with stepsize control
array LLG::BS23(const array& m, const double dt, double& err) {
    if (reject || calls == 0 || llg_wasnormalized) {
        heff = fheff(m);
        k[1] = fdmdt(m, heff);
    } else
        k[1] = k[4];

    heff = fheff(m + dt * (1. / 2. * k[1]));
    k[2] = fdmdt(m + dt * (1. / 2. * k[1]), heff);
    heff = fheff(m + dt * (+3. / 4. * k[2]));
    k[3] = fdmdt(m + dt * (+3. / 4. * k[2]), heff);
    sumbk = dt * (2. / 9. * k[1] + 1. / 3. * k[2] + 4. / 9. * k[3]);
    heff = fheff(m + sumbk);
    k[4] = fdmdt(m + sumbk, heff);

    rk_error = sumbk - dt * (7. / 24. * k[1] + 1. / 4. * k[2] + 1. / 3. * k[3] + 1. / 8. * k[4]);
    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));
    return sumbk;
}

// Standard Runge-Kutta 4th order method
array LLG::rk4(const array& m, const double dt) {
    heff = fheff(m);
    k1 = dt * fdmdt(m, heff);
    heff = fheff(m + 1. / 2. * k1);
    k2 = dt * fdmdt(m + 1. / 2. * k1, heff);
    heff = fheff(m + 1. / 2. * k2);
    k3 = dt * fdmdt(m + 1. / 2. * k2, heff);
    heff = fheff(m + k3);
    k4 = dt * fdmdt(m + k3, heff);

    return (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
}
// Standard Runge-Kutta 4th order method
array LLG::rk4minimal(const array& m, const double dt) {
    // array k, m_add;
    array k = array(state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    array m_add = array(state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);

    heff = fheff(m);
    k = dt * fdmdt(m, heff);
    m_add = k;
    heff = fheff(m + 1. / 2. * k);
    k = dt * fdmdt(m + 1. / 2. * k, heff);
    m_add += 2. * k;
    heff = fheff(m + 1. / 2. * k);
    k = dt * fdmdt(m + 1. / 2. * k, heff);
    m_add += 2. * k;
    heff = fheff(m + k);
    k = dt * fdmdt(m + k, heff);
    m_add += k;
    return m_add / 6.;
}

// Runge-Kutta 3/8 variation
array LLG::rk4_3o8(const array& m, const double dt) {
    heff = fheff(m);
    k1 = dt * fdmdt(m, heff);
    heff = fheff(m + 1. / 3. * k1);
    k2 = dt * fdmdt(m + 1. / 3. * k1, heff);
    heff = fheff(m - 1. / 3. * k1 + k2);
    k3 = dt * fdmdt(m - 1. / 3. * k1 + k2, heff);
    heff = fheff(m + k1 - k2 + k3);
    k4 = dt * fdmdt(m + k1 - k2 + k3, heff);

    return (k1 + 3. * k2 + 3. * k3 + k4) / 8.;
}
// Runge-Kutta-Fehlberg 5th order method
array LLG::RKF5(const array& m, const double dt) {

    heff = fheff(m);
    k1 = dt * fdmdt(m, heff);
    heff = fheff(m + 1. / 4. * k1);
    k2 = dt * fdmdt(m + 1. / 4. * k1, heff);
    heff = fheff(m + 3. / 32. * k1 + 9 / 32. * k2);
    k3 = dt * fdmdt(m + 3. / 32. * k1 + 9 / 32. * k2, heff);
    heff = fheff(m + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3);
    k4 = dt * fdmdt(m + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3, heff);
    heff = fheff(m + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4);
    k5 = dt * fdmdt(m + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4, heff);
    heff = fheff(m - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5);
    k6 = dt * fdmdt(m - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5, heff);

    sumbk = 16. / 135. * k1 + 6656. / 12825. * k3 + 28561. / 56430. * k4 - 9. / 50. * k5 + 2. / 55. * k6;
    return sumbk;
}

// Runge-Kutta-Fehlberg method with stepsize control
array LLG::RKF45(const array& m, const double dt, double& err) {
    // state0.t+=dt;

    heff = fheff(m);
    k1 = dt * fdmdt(m, heff);
    heff = fheff(m + 1. / 4. * k1);
    k2 = dt * fdmdt(m + 1. / 4. * k1, heff);
    heff = fheff(m + 3. / 32. * k1 + 9 / 32. * k2);
    k3 = dt * fdmdt(m + 3. / 32. * k1 + 9 / 32. * k2, heff);
    heff = fheff(m + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3);
    k4 = dt * fdmdt(m + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3, heff);
    heff = fheff(m + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4);
    k5 = dt * fdmdt(m + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4, heff);
    heff = fheff(m - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5);
    k6 = dt * fdmdt(m - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5, heff);

    sumbk = 16. / 135. * k1 + 6656. / 12825. * k3 + 28561. / 56430. * k4 - 9. / 50. * k5 + 2. / 55. * k6;
    rk_error = sumbk - (25. / 216. * k1 + 1408. / 2565. * k3 + 2197. / 4104. * k4 - 1. / 5. * k5);

    // rk_abs_error = maxnorm(rk_error);
    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));
    return sumbk;
}

// Cash-Karp  method
array LLG::CK45(const array& m, const double dt, double& err) {
    // state0.t+=dt;

    heff = fheff(m);
    k1 = dt * fdmdt(m, heff);

    rktemp = m + 1. / 5. * k1;
    heff = fheff(rktemp);
    k2 = dt * fdmdt(rktemp, heff);

    rktemp = m + 3. / 40. * k1 + 9. / 40. * k2;
    heff = fheff(rktemp);
    k3 = dt * fdmdt(rktemp, heff);

    rktemp = m + 3. / 10. * k1 - 9. / 10. * k2 + 6. / 5. * k3;
    heff = fheff(rktemp);
    k4 = dt * fdmdt(rktemp, heff);

    rktemp = m - 11. / 54. * k1 + 5. / 2. * k2 - 70. / 27. * k3 + 35. / 27. * k4;
    heff = fheff(rktemp);
    k5 = dt * fdmdt(rktemp, heff);

    rktemp =
        m + 1631. / 55296. * k1 + 175. / 512. * k2 + 575. / 13824. * k3 + 44275. / 110592. * k4 + 253. / 4096. * k5;
    heff = fheff(rktemp);
    k6 = dt * fdmdt(rktemp, heff);

    sumbk = 37. / 378. * k1 + 250. / 621. * k3 + 125. / 594. * k4 + 512. / 1771. * k6;
    rk_error =
        sumbk - (2825. / 27648. * k1 + 18575. / 48384. * k3 + 13525. / 55296. * k4 + 277. / 14336. * k5 + 1. / 4. * k6);

    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));
    // err = maxnorm(rk_error);
    return sumbk;
}

// Tsitorous 4/5
array LLG::tsit45(const array& m, const double dt, double& err) {
    if (reject || calls == 0 || llg_wasnormalized) {
        heff = fheff(m);
        k[1] = fdmdt(m, heff);
    } else
        k[1] = k[s];
    for (int i = 2; i <= s; i++) {
        rktemp = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
        for (int j = 1; j < i; j++) {
            rktemp += a[i][j] * k[j];
        }
        rktemp *= dt;
        heff = fheff(m + rktemp);
        k[i] = fdmdt(m + rktemp, heff);
    }
    // Local extrapolation using 5th order approx
    sumbk = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i < s; i++) {
        sumbk += a[s][i] * k[i];
    }
    sumbk *= dt;
    // Error estimation using 4th order approx
    rk_error = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i <= s; i++) {
        rk_error += bhat[i] * k[i];
    }
    rk_error *= dt;
    // rk_error=sumbk-rk_error;
    // print("rk_error", rk_error);
    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));
    return sumbk;
}

// Dormand-Prince 4/5 method
array LLG::DP45(const array& m, const double dt, double& err) {
    // Iterating over a-matrix and calculating k[i]s
    if (reject || calls == 0 || llg_wasnormalized) {
        heff = fheff(m);
        k[1] = fdmdt(m, heff);
    } else
        k[1] = k[s];
    for (int i = 2; i <= s; i++) {
        rktemp = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
        for (int j = 1; j < i; j++) {
            rktemp += a[i][j] * k[j];
        }
        rktemp *= dt;
        heff = fheff(m + rktemp);
        k[i] = fdmdt(m + rktemp, heff);
    }
    // Local extrapolation using 5th order approx
    sumbk = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i < s; i++) {
        sumbk += a[s][i] * k[i];
    }
    sumbk *= dt;
    // Error estimation using 4th order approx
    rk_error = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i <= s; i++) {
        rk_error += e[i] * k[i];
    }
    rk_error *= dt;
    //!!!Note: here e is already the difference between the ususal b and
    //! bhat!!!! (no rk_error=sumbk-rk_error)
    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));
    return sumbk;
}

// Bogacki 4, 5 method with sigle error andstepsize control
array LLG::BS45(const array& m, const double dt, double& err) {
    if (reject || calls == 0 || llg_wasnormalized) {
        // Note: in generalized rkcall use: if(FSAL==false || reject || calls==0
        // || llg_wasnormalized){
        heff = fheff(m);
        k[1] = fdmdt(m, heff);
    } else
        k[1] = k[s];
    for (int i = 2; i <= s; i++) {
        rktemp = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
        for (int j = 1; j < i; j++) {
            rktemp += a[i][j] * k[j];
        }
        rktemp *= dt;
        heff = fheff(m + rktemp);
        k[i] = fdmdt(m + rktemp, heff);
    }

    sumbk = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i < s; i++) {
        sumbk += a[s][i] * k[i];
    }
    sumbk *= dt;

    rk_error = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i <= s; i++) {
        rk_error += b[i] * k[i];
    }
    rk_error *= dt;
    rk_error = sumbk - rk_error;
    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));

    return sumbk;
}

// Bogacki 4, 5 method with double error estimation and with stepsize control
array LLG::BS45de(const array& m, const double dt, double& err) {
    double e[9];

    // C  The coefficients E(*) refer to an estimate of the local error based on
    // C  the first formula of order 4.  It is the difference of the fifth order
    // C  result, here located in A(8, *), and the fourth order result.  By
    // C  construction both E(2) and E(7) are zero.
    // C
    e[1] = -3.e0 / 1280.e0;
    e[2] = 0.e0;
    e[3] = 6561.e0 / 632320.e0;
    e[4] = -343.e0 / 20800.e0;
    e[5] = 243.e0 / 12800.e0;
    e[6] = -1.e0 / 95.e0;
    e[7] = 0.e0;

    if (reject || calls == 0 || llg_wasnormalized) {
        // Note: in generalized rkcall use: if(FSAL==false || reject || calls==0
        // || llg_wasnormalized){
        heff = fheff(m);
        k[1] = fdmdt(m, heff);
    } else
        k[1] = k[s];
    // Compute stages 2 to 6  and check error, if accepted, succeed to stage 8,
    // calc sumbk and check second error sumbk-bi*ki
    for (int i = 2; i <= 6; i++) { // i=2, ..., 6
        // for(int i=2;i<=s;i++){
        rktemp = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
        for (int j = 1; j < i; j++) {
            rktemp += a[i][j] * k[j];
        }
        rktemp *= dt;
        heff = fheff(m + rktemp);
        k[i] = fdmdt(m + rktemp, heff);
    }

    // Check first 4th order approx yielding directly the error
    rk_error = array(state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    rk_error = e[1] * k[1];
    for (int i = 3; i <= 6; i++) { // i=3, 4, 5, 6
        rk_error += e[i] * k[i];
    }
    rk_error *= dt;
    err = maxnorm(rk_error / controller.givescale(m)); // TODO we only use m, not max(m, m+sumbk)
                                                       // for error estimate! Is this convenient?
    // std::cout<<"Test"<<afvalue(rk_error(0, 0, 0,
    // 0))<<"\t"<<maxnorm(rk_error)<<"\t"<<err<<std::endl;
    if (err > 1.)
        std::cout << "RKErr>1" << std::endl;

    // We now check this error with the controller, if it passes (cont returns
    // true), hnext is changed temporary but then overwritten in step after
    // passing second controlling This 2 error method shoud save some
    // computational cost if the first error estimate already detects a too
    // large error, we save comp cost of 2 stages
    if (controller.success(err, h) == false) {
        std::cout << "CONTROOOL" << std::endl;
        return constant(10000000.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    }

    // 7th and 8th stage
    for (int i = 7; i <= 8; i++) {
        rktemp = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
        for (int j = 1; j < i; j++) {
            rktemp += a[i][j] * k[j];
        }
        rktemp *= dt;
        heff = fheff(m + rktemp);
        k[i] = fdmdt(m + rktemp, heff);
    }

    // Calc sumbk
    sumbk = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i < s; i++) {
        sumbk += a[s][i] * k[i];
    }
    sumbk *= dt;
    // Calc second and more precise 4th order error estimate
    rk_error = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i <= s; i++) {
        rk_error += b[i] * k[i];
    }
    rk_error *= dt;
    rk_error = sumbk - rk_error;
    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));

    return sumbk;
}

// Dormand-Prince 7/8 method with stepsize control
array LLG::DP78(const array& m, const double dt, double& err) {
    heff = fheff(m);
    k[1] = fdmdt(m, heff);
    // Iterating over a-matrix, calculating k[i]s
    for (int i = 2; i <= s; i++) {
        rktemp = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
        for (int j = 1; j < i; j++) {
            rktemp += a[i][j] * k[j];
        }
        rktemp *= dt;
        heff = fheff(m + rktemp);
        k[i] = fdmdt(m + rktemp, heff);
    }
    // Calculating 8th order approx (local extrapolation)
    sumbk = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i <= s; i++) {
        sumbk += bhat[i] * k[i];
    }
    sumbk *= dt;

    // Calculating 7th order approx for error approx
    rk_error = constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
    for (int i = 1; i <= s; i++) {
        rk_error += b[i] * k[i];
    }
    rk_error *= dt;
    rk_error = sumbk - rk_error;
    err = maxnorm(rk_error / controller.givescale(max(m, m + sumbk)));
    return sumbk;
}

LLG::LLG(State state0_in, std::vector<std::unique_ptr<FieldTerm>> Fieldterms_in)
    : Fieldterms(Fieldterms_in), state0(state0_in) {
    // LLG::LLG (State state0_in, double atol_in, double  rtol_in, double
    // hmax_in, double  hmin_in, std::vector<std::unique_ptr<FieldTerm> >
    // Fieldterms_in) : Fieldterms(Fieldterms_in), state0(state0_in),
    // atol(atol_in), rtol(rtol_in), hmax(hmax_in), hmin(hmin_in){

    if ((state0.material.mode == 0) || (state0.material.mode == 1) || (state0.material.mode == 2) ||
        (state0.material.mode == 3) || (state0.material.mode == 5)) {
        h = controller.hmax;
        hnext = controller.hmax;
    } else
        h = 1.01 * controller.hmin;

    // Tsit45
    if (state0.material.mode == 8) {
        s = 7;
        // const int s=7;
        // Declare arrays and init with zeros
        // double a[s+1][s+1], c[s+1], b[s+1], bhat[s+1];
        // for(int i=0;i<=s;i++){
        //  c[i]=0.;
        //  b[i]=0.;
        //  bhat[i]=0.;
        //  for(int j=0;j<=s;j++){
        //    a[i][j]=0.;
        //  }
        //}
        c[2] = 0.161;
        c[4] = 0.9;
        c[6] = 1.;
        c[7] = 1.;
        b[2] = 0.01;
        b[4] = 1.379008574103742;
        b[6] = 2.324710524099774;
        bhat[2] = 0.000816434459657;
        bhat[4] = 0.144711007173263;
        bhat[6] = 0.458082105929187;
        a[3][2] = 0.3354806554923570;
        a[5][2] = -11.74888356406283;
        a[5][3] = 7.495539342889836;
        a[6][2] = -12.92096931784711;
        a[6][4] = -0.07158497328140100;
        c[3] = 0.327;
        c[5] = 0.9800255409045097;
        b[1] = 0.09646076681806523;
        b[3] = 0.4798896504144996;
        b[5] = -3.290069515436081;
        bhat[1] = 0.001780011052226;
        bhat[3] = -0.007880878010262;
        bhat[5] = -0.582357165452555;
        bhat[7] = 1. / 66.;
        b[7] = 0;
        a[4][2] = -6.359448489975075;
        a[4][3] = 4.362295432869581;
        a[5][4] = -0.09249506636175525;
        a[6][3] = 8.159367898576159;
        a[6][5] = -0.02826905039406838;
        a[2][1] = c[2];
        a[3][1] = c[3] - a[3][2];
        a[4][1] = c[4] - a[4][2] - a[4][3];
        a[5][1] = c[5] - a[5][2] - a[5][3] - a[5][4];
        a[6][1] = c[6] - a[6][2] - a[6][3] - a[6][4] - a[6][5];
        for (int i = 1; i <= 6; i++) {
            a[7][i] = b[i];
        }
    }
    // DP45
    if (state0.material.mode == 9) {
        s = 7;
        // Declare arrays and init with zeros
        // double a[s+1][s+1], e[s+1];//, c[s+1]
        // for(int i=0;i<=s;i++){
        //  //c[i]=0.;
        //  e[i]=0.;
        //  for(int j=0;j<=s;j++){
        //    a[i][j]=0.;
        //  }
        //}

        // c[2]=0.2, c[3]=0.3, c[4]=0.8, c[5]=8.0/9.0;
        a[2][1] = 0.2, a[3][1] = 3.0 / 40.0, a[3][2] = 9.0 / 40.0, a[4][1] = 44.0 / 45.0, a[4][2] = -56.0 / 15.0,
        a[4][3] = 32.0 / 9.0, a[5][1] = 19372.0 / 6561.0, a[5][2] = -25360.0 / 2187.0, a[5][3] = 64448.0 / 6561.0,
        a[5][4] = -212.0 / 729.0, a[6][1] = 9017.0 / 3168.0, a[6][2] = -355.0 / 33.0, a[6][3] = 46732.0 / 5247.0,
        a[6][4] = 49.0 / 176.0, a[6][5] = -5103.0 / 18656.0, a[7][1] = 35.0 / 384.0, a[7][3] = 500.0 / 1113.0,
        a[7][4] = 125.0 / 192.0, a[7][5] = -2187.0 / 6784.0, a[7][6] = 11.0 / 84.0, e[1] = 71.0 / 57600.0,
        e[3] = -71.0 / 16695.0, e[4] = 71.0 / 1920.0, e[5] = -17253.0 / 339200.0, e[6] = 22.0 / 525.0,
        e[7] = -1.0 / 40.0;
    }
    // BS5 and BS5de
    if (state0.material.mode == 10 || state0.material.mode == 11) {
        s = 8;
        FSAL = true;
        a[2][1] = 1.0e0 / 6.0e0;
        a[3][1] = 2.e0 / 27.e0;
        a[3][2] = 4.e0 / 27.e0;
        a[4][1] = 183.e0 / 1372.e0;
        a[4][2] = -162.e0 / 343.e0;
        a[4][3] = 1053.e0 / 1372.e0;
        a[5][1] = 68.e0 / 297.e0;
        a[5][2] = -4.e0 / 11.e0;
        a[5][3] = 42.e0 / 143.e0;
        a[5][4] = 1960.e0 / 3861.e0;
        a[6][1] = 597.e0 / 22528.e0;
        a[6][2] = 81.e0 / 352.e0;
        a[6][3] = 63099.e0 / 585728.e0;
        a[6][4] = 58653.e0 / 366080.e0;
        a[6][5] = 4617.e0 / 20480.e0;
        a[7][1] = 174197.e0 / 959244.e0;
        a[7][2] = -30942.e0 / 79937.e0;
        a[7][3] = 8152137.e0 / 19744439.e0;
        a[7][4] = 666106.e0 / 1039181.e0;
        a[7][5] = -29421.e0 / 29068.e0;
        a[7][6] = 482048.e0 / 414219.e0;
        a[8][1] = 587.e0 / 8064.e0;
        a[8][2] = 0.e0;
        a[8][3] = 4440339.e0 / 15491840.e0;
        a[8][4] = 24353.e0 / 124800.e0;
        a[8][5] = 387.e0 / 44800.e0;
        a[8][6] = 2152.e0 / 5985.e0;
        a[8][7] = 7267.e0 / 94080.e0;
        // C
        // C
        c[1] = 0.e0;
        c[2] = 1.e0 / 6.e0;
        c[3] = 2.e0 / 9.e0;
        c[4] = 3.e0 / 7.e0;
        c[5] = 2.e0 / 3.e0;
        c[6] = 3.e0 / 4.e0;
        c[7] = 1.e0;
        c[8] = 1.e0;
        // C  The coefficients B(*) refer to the formula of order 4.
        b[1] = 2479.e0 / 34992.e0;
        b[2] = 0.e0;
        b[3] = 123.e0 / 416.e0;
        b[4] = 612941.e0 / 3411720.e0;
        b[5] = 43.e0 / 1440.e0;
        b[6] = 2272.e0 / 6561.e0;
        b[7] = 79937.e0 / 1113912.e0;
        b[8] = 3293.e0 / 556956.e0;
    }
    if (state0.material.mode == 12) {
        s = 13;
        FSAL = false;
        a[2][1] = 5.55555555555555555555555555556e-2;
        a[3][1] = 2.08333333333333333333333333333e-2;
        a[3][2] = 6.25e-2;
        a[4][1] = 3.125e-2;
        a[4][2] = 0.e0;
        a[4][3] = 9.375e-2;
        a[5][1] = 3.125e-1;
        a[5][2] = 0.e0;
        a[5][3] = -1.171875e0;
        a[5][4] = 1.171875e0;
        a[6][1] = 3.75e-2;
        a[6][2] = 0.e0;
        a[6][3] = 0.e0;
        a[6][4] = 1.875e-1;
        a[6][5] = 1.5e-1;
        a[7][1] = 4.79101371111111111111111111111e-2;
        a[7][2] = 0.e0;
        a[7][3] = 0.0e0;
        a[7][4] = 1.12248712777777777777777777778e-1;
        a[7][5] = -2.55056737777777777777777777778e-2;
        a[7][6] = 1.28468238888888888888888888889e-2;
        a[8][1] = 1.6917989787292281181431107136e-2;
        a[8][2] = 0.e0;
        a[8][3] = 0.e0;
        a[8][4] = 3.87848278486043169526545744159e-1;
        a[8][5] = 3.59773698515003278967008896348e-2;
        a[8][6] = 1.96970214215666060156715256072e-1;
        a[8][7] = -1.72713852340501838761392997002e-1;
        a[9][1] = 6.90957533591923006485645489846e-2;
        a[9][2] = 0.e0;
        a[9][3] = 0.e0;
        a[9][4] = -6.34247976728854151882807874972e-1;
        a[9][5] = -1.61197575224604080366876923982e-1;
        a[9][6] = 1.38650309458825255419866950133e-1;
        a[9][7] = 9.4092861403575626972423968413e-1;
        a[9][8] = 2.11636326481943981855372117132e-1;
        a[10][1] = 1.83556996839045385489806023537e-1;
        a[10][2] = 0.e0;
        a[10][3] = 0.e0;
        a[10][4] = -2.46876808431559245274431575997e0;
        a[10][5] = -2.91286887816300456388002572804e-1;
        a[10][6] = -2.6473020233117375688439799466e-2;
        a[10][7] = 2.84783876419280044916451825422e0;
        a[10][8] = 2.81387331469849792539403641827e-1;
        a[10][9] = 1.23744899863314657627030212664e-1;
        a[11][1] = -1.21542481739588805916051052503e0;
        a[11][2] = 0.e0;
        a[11][3] = 0.e0;
        a[11][4] = 1.66726086659457724322804132886e1;
        a[11][5] = 9.15741828416817960595718650451e-1;
        a[11][6] = -6.05660580435747094755450554309e0;
        a[11][7] = -1.60035735941561781118417064101e1;
        a[11][8] = 1.4849303086297662557545391898e1;
        a[11][9] = -1.33715757352898493182930413962e1;
        a[11][10] = 5.13418264817963793317325361166e0;
        a[12][1] = 2.58860916438264283815730932232e-1;
        a[12][2] = 0.e0;
        a[12][3] = 0.e0;
        a[12][4] = -4.77448578548920511231011750971e0;
        a[12][5] = -4.3509301377703250944070041181e-1;
        a[12][6] = -3.04948333207224150956051286631e0;
        a[12][7] = 5.57792003993609911742367663447e0;
        a[12][8] = 6.15583158986104009733868912669e0;
        a[12][9] = -5.06210458673693837007740643391e0;
        a[12][10] = 2.19392617318067906127491429047e0;
        a[12][11] = 1.34627998659334941535726237887e-1;
        a[13][1] = 8.22427599626507477963168204773e-1;
        a[13][2] = 0.e0;
        a[13][3] = 0.e0;
        a[13][4] = -1.16586732572776642839765530355e1;
        a[13][5] = -7.57622116690936195881116154088e-1;
        a[13][6] = 7.13973588159581527978269282765e-1;
        a[13][7] = 1.20757749868900567395661704486e1;
        a[13][8] = -2.12765911392040265639082085897e0;
        a[13][9] = 1.99016620704895541832807169835e0;
        a[13][10] = -2.34286471544040292660294691857e-1;
        a[13][11] = 1.7589857770794226507310510589e-1;
        a[13][12] = 0.e0;
        // C
        // C  The coefficients BHAT(*) refer to the formula used to advance the
        // C  integration, here the one of order 8.  The coefficients B(*) refer
        // C  to the other formula, here the one of order 7.
        // C
        bhat[1] = 4.17474911415302462220859284685e-2;
        bhat[2] = 0.e0;
        bhat[3] = 0.e0;
        bhat[4] = 0.e0;
        bhat[5] = 0.e0;
        bhat[6] = -5.54523286112393089615218946547e-2;
        bhat[7] = 2.39312807201180097046747354249e-1;
        bhat[8] = 7.0351066940344302305804641089e-1;
        bhat[9] = -7.59759613814460929884487677085e-1;
        bhat[10] = 6.60563030922286341461378594838e-1;
        bhat[11] = 1.58187482510123335529614838601e-1;
        bhat[12] = -2.38109538752862804471863555306e-1;
        bhat[13] = 2.5e-1;
        // C
        b[1] = 2.9553213676353496981964883112e-2;
        b[2] = 0.e0;
        b[3] = 0.e0;
        b[4] = 0.e0;
        b[5] = 0.e0;
        b[6] = -8.28606276487797039766805612689e-1;
        b[7] = 3.11240900051118327929913751627e-1;
        b[8] = 2.46734519059988698196468570407e0;
        b[9] = -2.54694165184190873912738007542e0;
        b[10] = 1.44354858367677524030187495069e0;
        b[11] = 7.94155958811272872713019541622e-2;
        b[12] = 4.44444444444444444444444444445e-2;
        b[13] = 0.e0;
        // C
        c[1] = 0.e0;
        c[2] = 5.55555555555555555555555555556e-2;
        c[3] = 8.33333333333333333333333333334e-2;
        c[4] = 1.25e-1;
        c[5] = 3.125e-1;
        c[6] = 3.75e-1;
        c[7] = 1.475e-1;
        c[8] = 4.65e-1;
        c[9] = 5.64865451382259575398358501426e-1;
        c[10] = 6.5e-1;
        c[11] = 9.24656277640504446745013574318e-1;
        c[12] = 1.e0;
        c[13] = c[12];
    }
}

// array LLG::llgstepEEold(array& m, double dt){
//    timer_integrator = timer::start();
//    heff = Demag.solve(m) + Exch.solve(m);
//    crosstemp  =  cross4(m, heff);
//    array dmdt = - state0.material.gamma/(1.+pow(state0.material.alpha, 2)) *
//    cross4(m, heff) -
//    state0.material.alpha*state0.material.gamma/(1.+pow(state0.material.alpha,
//    2)) * cross4(m, crosstemp); m += dt * dmdt;
//    if(state0.state.material.afsync) af::sync();
//    time_integrator += timer::stop(timer_integrator);
////    std::cout << "T: heff.dims " << heff.dims() << " m.dims " << m.dims() <<
///" cross4 " << crosstemp.dims()  << std::endl;
//    return m/tile(sqrt(sum(m*m, 3)), 1, 1, 1, 3);
//};
// array LLG::llgstepEEold(array& m, double dt, array& h_zee){
//    timer_integrator = timer::start();
//    heff = Demag.solve(m) + Exch.solve(m) + h_zee;
//    crosstemp  =  cross4(m, heff); array dmdt = -
//    state0.material.gamma/(1.+pow(state0.material.alpha, 2)) * cross4(m, heff)
//    -
//    state0.material.alpha*state0.material.gamma/(1.+pow(state0.material.alpha,
//    2)) * cross4(m, crosstemp); m += dt * dmdt;
//    if(state0.state.material.afsync) af::sync(); time_integrator +=
//    timer::stop(timer_integrator); return m/tile(sqrt(sum(m*m, 3)), 1, 1, 1,
//    3); };
//
////Explicit Euler LLG stpe without zeeman field
// array LLG::llgstepEE(array& m){
//    timer_integrator = timer::start();
//    m += explicitEuler(m, state0.material.dt);
//    if(state0.state.material.afsync) af::sync();
//    time_integrator += timer::stop(timer_integrator);
//    return m/tile(sqrt(sum(m*m, 3)), 1, 1, 1, 3);
//};

//// Dormand-Prince 4/5 method
// array LLG::DP45(const array& m, const double dt, double& err)
//{
//  //double c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;
//  double a21=0.2, a31=3.0/40.0, a32=9.0/40.0, a41=44.0/45.0, a42=-56.0/15.0,
//  a43=32.0/9.0, a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0,
//  a54=-212.0/729.0, a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0,
//  a64=49.0/176.0, a65=-5103.0/18656.0, a71=35.0/384.0, a73=500.0/1113.0,
//  a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0, e1=71.0/57600.0,
//  e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0,
//  e7=-1.0/40.0; const double e[8]={0, e1, 0, e3, e4, e5, e6, e7}; const double
//  a[8][8]={
//    {0, 0, 0, 0, 0, 0, 0, 0},
//    {0, 0, 0, 0, 0, 0, 0, 0},
//    {0, a21, 0, 0, 0, 0, 0, 0},
//    {0, a31, a32, 0, 0, 0, 0, 0},
//    {0, a41, a42, a43, 0, 0, 0, 0},
//    {0, a51, a52, a53, a54, 0, 0, 0},
//    {0, a61, a62, a63, a64, a65, 0, 0},
//    {0, a71, 0  , a73, a74, a75, a76, 0.}
//  };
//
// // if(reject) std::cout<<"!!!!!!!! Prev was rejected"<<std::endl;
////  std::cout<<"mini = "<<afvalue((m)(0, 0, 0, 0))<<std::endl;
//
//  if(reject || calls==0 || llg_wasnormalized){
//    heff =  fheff(m);
//    k1   =  fdmdt(m, heff);
//    }
//  else
//    k1=k7;
//
////  std::cout<<"h of k1 = "<<afvalue((heff)(0, 0, 0, 0))<<"\n"<<std::endl;
////  std::cout<<"k1 = "<<afvalue((k1)(0, 0, 0, 0))<<"\n"<<std::endl;
////
////  std::cout<<"m k1 = "<<afvalue((m)(0, 0, 0, 0))<<"\t"<<"h of k1 =
///"<<afvalue((heff)(0, 0, 0, 0))<<"\t"<<"k1= "<<afvalue(k1(0, 0, 0, 0))<<"\t
/// dtfdmdt= "<<afvalue(fdmdt(m, heff)(0, 0, 0, 0))<<"\n"<<std::endl;
//  heff =  fheff(m + dt *( a[2][1] * k1) ); k2   =  fdmdt(m + dt *( a[2][1] *
//  k1) , heff); heff =  fheff(m + dt *( a[3][1] * k1 + a[3][2] * k2) ); k3   =
//  fdmdt(m + dt *( a[3][1] * k1 + a[3][2] * k2) , heff); heff =  fheff(m + dt
//  *( a[4][1] * k1 + a[4][2] * k2 +  a[4][3] * k3) ); k4   =  fdmdt(m + dt *(
//  a[4][1] * k1 + a[4][2] * k2 +  a[4][3] * k3) , heff); heff =  fheff(m + dt
//  *( a[5][1] * k1 + a[5][2] * k2 +  a[5][3] * k3 + a[5][4] * k4) ); k5   =
//  fdmdt(m + dt *( a[5][1] * k1 + a[5][2] * k2 +  a[5][3] * k3 + a[5][4] * k4)
//  , heff); heff =  fheff(m + dt *( a[6][1] * k1 + a[6][2] * k2 +  a[6][3] * k3
//  + a[6][4] * k4 + a[6][5] * k5)                     ); k6   =  fdmdt(m + dt
//  *( a[6][1] * k1 + a[6][2] * k2 +  a[6][3] * k3 + a[6][4] * k4 + a[6][5] *
//  k5)              , heff); heff =  fheff(m + dt *( a[7][1] * k1 +  a[7][3] *
//  k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6)      ); k7   =  fdmdt(m +
//  dt *( a[7][1] * k1                +  a[7][3] * k3 + a[7][4] * k4 + a[7][5] *
//  k5 + a[7][6] * k6), heff);
//
//  sumbk    = dt*( a[7][1]*k1 + a[7][3]*k3 + a[7][4]*k4 + a[7][5]*k5 +
//  a[7][6]*k6); rk_error = sumbk - dt*(e[1]*k1 + e[2]*k2 + e[3]*k3 + e[4]*k4 +
//  e[5]*k5 + e[6]*k6 + e[7]*k7); err=maxnorm(rk_error/givescale(max(m,
//  m+sumbk)));
//  //std::cout<<"mk7  = "<<afvalue((m + a[7][1] * k1                +  a[7][3]
//  * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6)(0, 0, 0, 0))<<std::endl;
////  std::cout<<"after heff = "<<afvalue((m + a[7][1] * k1                +
/// a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6 )(0, 0, 0,
/// 0))<<std::endl; /  std::cout<<"after k7   = "<<afvalue((m + a[7][1] * k1 +
/// a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6 )(0, 0, 0,
/// 0))<<"\n"<<std::endl; /  std::cout<<"h of k7 = "<<afvalue((heff)(0, 0, 0,
/// 0))<<std::endl;
//
//  //Todo: not differs in 7th digit
//  //std::cout.precision(12);
//  //std::cout << maxnorm(k1(0, 0, 0, 0)) << "\t" << maxnorm(k7(0, 0, 0, 0)) <<
//  std::endl;
////  std::cout<<"m+sumbk = "<<afvalue((m+sumbk)(0, 0, 0, 0))<<std::endl;
//  //rk_error = (a[7][1] - e[1] ) * k1 + (a[7][2] - e[2] * k2 ) + (a[7][3] -
//  e[3] ) * k3 + ( a[7][4] -e[4] ) * k4 + (a[7][5] - e[5]) * k5 + (a[7][6] -
//  e[6] )* k6 -e[7] * k7;
//  //std::cout << "h = "<<h<< " error = "<<err<<" maxnorm(rk_error)
//  "<<maxnorm(rk_error)<<" givescale "<< maxnorm(givescale(max(m,
//  m+sumbk)))<<std::endl;
//  //rk_abs_error = maxnorm(rk_error);
//  //std::cout<<"m+s3 = "<<afvalue((m+sumbk)(0, 0, 0, 0))<<std::endl;
//
////  std::cout<<"k7 = "<<afvalue((dt * fdmdt(m + a[7][1] * k1                +
/// a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6, heff))(0, 0, 0,
/// 0))<<std::endl;
//  //std::cout<<"m k7 = "<<afvalue((m+sumbk)(0, 0, 0, 0))<<"\t"<<"h of k7 =
//  "<<afvalue((heff)(0, 0, 0, 0))<<"\t"<<"k7= "<<afvalue(k7(0, 0, 0,
//  0))<<std::endl;
//
////  std::cout<<"m k7 = "<<afvalue((m + dt*(a[7][1] * k1+  a[7][3] * k3 +
/// a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6))(0, 0, 0, 0))<<"\t"<<"h of k7 =
///"<<afvalue((heff)(0, 0, 0, 0))<<"\t"<<"k7= "<<afvalue(k7(0, 0, 0, 0))<<"\t
/// dtfdmdt= "<<afvalue(fdmdt(m + dt*(a[7][1] * k1                +  a[7][3] *
/// k3
///+ a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6), heff)(0, 0, 0, 0))<<std::endl;
//  return sumbk;
//}

//// Dormand-Prince 4/5 method
// array LLG::DP45(const array& m, const double dt, double& err)
//{
//  if(reject || calls==0 || llg_wasnormalized){
//    heff =  fheff(m);
//    k1   =  fdmdt(m, heff);
//    }
//  else
//    k1=k7;
//
//  //std::cout<<"m k1 = "<<afvalue((m)(0, 0, 0, 0))<<"\t"<<"h of k1 =
//  "<<afvalue((heff)(0, 0, 0, 0))<<"\t"<<"k1= "<<afvalue(k1(0, 0, 0, 0))<<"\t
//  dtfdmdt= "<<afvalue(fdmdt(m, heff)(0, 0, 0, 0))<<"\n"<<std::endl; rktemp = m
//  + dt * (1./5. *k1);
//
//  heff = fheff(rktemp      );
//  k2   = fdmdt(rktemp, heff);
//
//  rktemp = m + dt * ( 3./40. *k1 + 9./40. * k2);
//
//  heff = fheff(rktemp      );
//  k3   = fdmdt(rktemp, heff);
//
//  rktemp    = m + dt * (44./45. * k1 - 56./15. *k2 + 32./9. *k3);
//
//  heff = fheff(rktemp      );
//  k4   = fdmdt(rktemp, heff);
//
//  rktemp    = m + dt * (19372./6561. *k1 - 25360./2187. *k2 + 64448./6561. *
//  k3 - 212./729. * k4);
//
//  heff = fheff(rktemp      );
//  k5   = fdmdt(rktemp, heff);
//
//  rktemp    = m + dt * (9017./3168. * k1 -355./33. * k2 + 46732./5247. * k3
//  + 49./176. * k4 -5103./18656. * k5);
//
//  heff = fheff(rktemp      );
//  k6   = fdmdt(rktemp, heff);
//
//  rktemp    = m + dt * (35./384. * k1                     + 500./1113. * k3 +
//  125./192. * k4 -2187./6784. * k5 + 11./84. * k6);
//
//  heff = fheff(rktemp      );
//  k7   = fdmdt(rktemp, heff);
//
//  sumbk =        dt * ( 35./384. * k1                    + 500./1113. * k3 +
//  125./192. * k4 -2187./6784. * k5 + 11./84. * k6); rk_error = sumbk - dt * (
//  5179./57600. * k1           + 7571./16695. * k3 + 393/640. * k4
//  -92097./339200.* k5 + 187./2100. * k6 + 1./40. * k7);
//  //std::cout<<"m k7 = "<<afvalue((m)(0, 0, 0, 0))<<"\t"<<"h of k7 =
//  "<<afvalue((heff)(0, 0, 0, 0))<<"\t"<<"k7= "<<afvalue(k7(0, 0, 0, 0))<<"\t
//  dtfdmdt= "<<afvalue(fdmdt(m+sumbk, heff)(0, 0, 0, 0))<<std::endl;
//
//  err=maxnorm(rk_error/givescale(max(m, m+sumbk)));
//  return sumbk;
//}

//// Bogacki 4, 5 method with stepsize control
// array LLG::BS45(const array& m, const double h , double& err)
//{
//
//  heff =       fheff(m ); k1   =  h  * fdmdt(m , heff); heff =       fheff(m
//  +  1.0e0/6.0e0         * k1 ); k2   =  h  * fdmdt(m   +  1.0e0/6.0e0 * k1 ,
//  heff); heff =       fheff(m   +  2.e0/27.e0          * k1  +  4.e0/27.e0
//  * k2 ); k3   =  h  * fdmdt(m   +  2.e0/27.e0          * k1  +  4.e0/27.e0
//  * k2 , heff); heff =       fheff(m   +  183.e0/1372.e0      * k1
//  -162.e0/343.e0    * k2   + 1053.e0/1372.e0       * k3 ); k4   =  h  *
//  fdmdt(m   +  183.e0/1372.e0      * k1  -162.e0/343.e0    * k2   +
//  1053.e0/1372.e0       * k3 , heff); heff =       fheff(m   +  68.e0/297.e0
//  * k1  -4.e0/11.e0       * k2   + 42.e0/143.e0          * k3  +
//  1960.e0/3861.e0     * k4 ); k5   =  h  * fdmdt(m   +  68.e0/297.e0        *
//  k1  -4.e0/11.e0       * k2   + 42.e0/143.e0          * k3  + 1960.e0/3861.e0
//  * k4                                                                       ,
//  heff); heff =       fheff(m   +  597.e0/22528.e0     * k1  +  81.e0/352.e0
//  * k2   +  63099.e0/585728.e0   * k3  + 58653.e0/366080.e0  * k4 +
//  4617.e0/20480.e0 * k5                                                     );
//  k6   =  h  * fdmdt(m   +  597.e0/22528.e0     * k1  +  81.e0/352.e0   * k2
//  +  63099.e0/585728.e0   * k3  + 58653.e0/366080.e0  * k4 + 4617.e0/20480.e0
//  * k5                                               , heff); heff = fheff(m
//  +  174197.e0/959244.e0 * k1  -30942.e0/79937.e0* k2 +8152137.e0/19744439.e0
//  * k3  + 666106.e0/1039181.e0* k4 - 29421.e0/29068.e0* k5 +
//  482048.e0/414219.e0*k6                            ); k7   =  h  * fdmdt(m +
//  174197.e0/959244.e0 * k1  -30942.e0/79937.e0* k2   +8152137.e0/19744439.e0 *
//  k3  + 666106.e0/1039181.e0* k4 - 29421.e0/29068.e0* k5 +
//  482048.e0/414219.e0*k6                      , heff); heff =       fheff(m +
//  587.e0/8064.e0      * k1                           + 4440339.e0/15491840.e0*
//  k3  + 24353.e0/124800.e0  * k4 + 387.e0/44800.e0 * k5 + 2152.e0/5985.e0 *k6
//  + 7267.e0/94080.e0*k7      ); k8   =  h  * fdmdt(m   +  587.e0/8064.e0 * k1
//  + 4440339.e0/15491840.e0* k3  + 24353.e0/124800.e0  * k4 + 387.e0/44800.e0 *
//  k5 + 2152.e0/5985.e0    *k6 + 7267.e0/94080.e0*k7, heff);
//
//  //c instead of b: sumbk     =        1.e0/6.e0 *k2  +2.e0/9.e0 * k3
//  +3.e0/7.e0 * k4 +2.e0/3.e0   * k5  + 3.e0/4.e0*k6 + 1.e0 *k7   +1.e0 * k8 ;
//  sumbk     =       2479.e0/34992.e0 * k1  + 123.e0/416.e0 * k3  +
//  612941.e0/3411720.e0 * k4 +43.e0/1440.e0 * k5  + 2272.e0/6561.e0*k6 +
//  79937.e0/1113912.e0 *k7   +3293.e0/556956.e0* k8    ; rk_error = sumbk - (
//  -3.e0/1280.e0  * k1  + 6561.e0/632320.e0* k3   -343.e0/20800.e0* k4
//  +243.e0/12800.e0  * k5  -1.e0/95.e0*k6);
//
//  err=maxnorm(rk_error/givescale(max(m, m+sumbk)));
//  //err = maxnorm(rk_error);
//  return sumbk;
//}

//// Runge-Kutta-Fehlberg method with stepsize control
// array LLG::RKF45(array m, double dt, double& rk_abs_error)
//{
//  //state0.t+=dt;
//
//  heff =       fheff(m ); k1   =  dt * fdmdt(m , heff);
//
//  k2   =  m   +    1./4.    * k1 ; heff =      fheff(k2       ); k2   = dt *
//  fdmdt(k2 , heff);
//
//  k3   =  m   +    3./32.   * k1  + 9/32.       * k2 ; heff =      fheff(k3 );
//  k3   = dt * fdmdt(k3               , heff);
//
//  k4   =  m   + 1932./2197. * k1  - 7200./2197. * k2   +  7296./2197. * k3 ;
//  heff =      fheff(k4                     );
//  k4   = dt * fdmdt(k4               , heff);
//
//  k5   =  m   +  439./216.  * k1  -     8.      * k2   +  3680./513.  * k3  -
//  845./4104. * k4                               ; heff =      fheff(k5 ); k5
//  = dt * fdmdt(k5               , heff);
//
//  k6   =  m   -    8./27.   * k1  +     2.      * k2   -  3544./2565. * k3  +
//  1859./4104. * k4  - 11./40. * k5               ; heff =      fheff(k6 ); k6
//  = dt * fdmdt(k6               , heff);
//
//  sumbk     =                16./135.  * k1  + 6656./12825.* k3  +
//  28561./56430.* k4    -9./50. * k5 + 2./55. *k6         ; rk_error = sumbk -
//  (      25./216.  * k1  +                       1408./2565. * k3  +
//  2197./4104. * k4    -1./5.  * k5                     );
//
//  double *rk_abs_error_host=NULL;
//  rk_abs_error_host = max(max(max(max(abs(rk_error), 0), 1), 2),
//  3).host<double>(); rk_abs_error = rk_abs_error_host[0];
//  freeHost(rk_abs_error_host);
//  return sumbk;
//}

//// Dormand-Prince 4/5 method
// array LLG::DP45(array m, double dt, double& rk_abs_error)
//{
//  //state0.t+=dt;
//  //if(calls == 0) {
//  heff =      fheff(m ); k1   = dt * fdmdt(m , heff);
//  //std::cout << "TEST: calls = 0" << std::endl;
//  //}
//  //else {
//  //  k1=k7;
//  //}
//
//  k2   = m + 1./5. *k1;
//  heff =      fheff(k2 ); k2   = dt * fdmdt(k2 , heff);
//
//  k3   = m + 3./40. *k1 + 9./40. * k2;
//  heff =      fheff(k3                     );
//  k3   = dt * fdmdt(k3               , heff);
//
//  k4   = m + 44./45. * k1 - 56./15. *k2 + 32./9. *k3;
//  heff =      fheff(k4                     );
//  k4   = dt * fdmdt(k4               , heff);
//
//  k5   = m + 19372./6561. *k1 - 25360./2187. *k2 + 64448./6561. * k3 -
//  212./729. * k4; heff =      fheff(k5                     ); k5   = dt *
//  fdmdt(k5               , heff);
//
//  k6   = m + 9017./3168. * k1 -355./33. * k2 + 46732./5247. * k3 + 49./176. *
//  k4 -5103./18656. * k5; heff =      fheff(k6                     ); k6   = dt
//  * fdmdt(k6               , heff);
//
//  k7   = m + 35./384. * k1 + 500./1113. * k3 + 125./192. * k4 -2187./6784. *
//  k5 + 11./84. * k6; heff =      fheff(k7                     ); k7   = dt *
//  fdmdt(k7               , heff);
//
//  sumbk =    35./384. * k1 + 500./1113. * k3 + 125./192. * k4 -2187./6784. *
//  k5 + 11./84. * k6; rk_error = sumbk - ( 5179./57600. * k1 + 7571./16695. *
//  k3 + 393/640. * k4 -92097./339200.* k5 + 187./2100. * k6 + 1./40. * k7);
//
//  double *rk_abs_error_host=NULL;
//  rk_abs_error_host = max(max(max(max(abs(rk_error), 0), 1), 2),
//  3).host<double>(); rk_abs_error = rk_abs_error_host[0];
//  freeHost(rk_abs_error_host);
//  return sumbk;
//}
//

////RK4
// array LLG::step(array& m){
//    timer_integrator = timer::start();
//    m += rk4(m, state0.material.dt);
//    if(state0.state.material.afsync) af::sync();
//    time_integrator += timer::stop(timer_integrator);
//    mean(m); //TODO this is needed to avoid cuda crash?!?
//    return m/tile(sqrt(sum(m*m, 3)), 1, 1, 1, 3);
//};

////Explicit Euler LLG stpe without zeeman field
// array LLG::llgstepEtesting(array& m, double dt){
//    timer_integrator = timer::start();
//    m += explicitEuler(m, dt);
//    if(state0.state.material.afsync) af::sync();
//    time_integrator += timer::stop(timer_integrator);
//    std::cout << "TEST:bef" << " m.dims = "  << m.dims() << " heff.dims = "<<
//    heff.dims() << std::endl; m/=tile(sqrt(sum(m*m, 3)), 1, 1, 1, 3);
//    std::cout << "TEST:aft" << " m.dims = "  << m.dims() << " heff.dims = "<<
//    heff.dims() << std::endl; return m;
//};

//    else if (state0.material.mode == 3){
//      double *mdiff=NULL;
//      array mtemp = m;
//      //m += RKF45(m, state0.material.dt, rk_abs_error);
//      m += RKF45(m, h, rk_abs_error);
//      mdiff=max(max(max(max(abs(m-mtemp), 0), 1), 2), 3).host<double>();
//      //mdiff=max(abs(m-mtemp)).host<double>();
//      rk_rel_error=rk_abs_error/(mdiff[0]);
//      freeHost(mdiff);
//      //For fixed stepsize: h_abs = state0.material.dt *
//      pow(headroom*rk_abs_tol_error/rk_abs_error, (1./(6.+1.)));
//      //For fixed stepsize: h_rel = state0.material.dt *
//      pow(headroom*rk_rel_tol_error/rk_rel_error, (1./(6.))); h_abs = h *
//      pow(headroom*rk_abs_tol_error/rk_abs_error, (1./(6.+1.))); //TODO is
//      s=6? h_rel = h * pow(headroom*rk_rel_tol_error/rk_rel_error, (1./(6.)));
//      //TODO is s=6? (h_abs < h_rel) ? h=h_abs : h=h_rel ;
//      std::cout.precision(12);
//      std::cout << "rk_abs_error : "<< rk_abs_error << "\t" << " rk_rel_error:
//      " << rk_rel_error  << "\t"
//       <<"h_abs "<<  h_abs  << "\t"<<  " h_rel= " << h_rel  << "\t"<< " h= "
//       << h  <<std::endl;
//    }

// array LLG::rk4s(array f(array, array), double dt, array m, array heff)
////array LLG::rk4s(array (*f)(array, array), double dt, array m, array heff)
//{
//	array	k1 = dt * f(m,          heff),
//		k2 = dt * f(m + k1 / 2, heff),
//		k3 = dt * f(m + k2 / 2, heff),
//		k4 = dt * f(m + k3,     heff);
//	return  (k1 + 2 * k2 + 2 * k3 + k4) / 6;
//}

// Slower: factor 1.14 :relative 1.2 instead of 1.06
// array LLG::cross4(array& a, array& b){
//  a=reorder(a, 3, 0, 1, 2);
//  //array a=reorder(ain, 3, 0, 1, 2);
//  b=reorder(b, 3, 0, 1, 2);
//  array c=constant(0.0, 3, a.dims(1), a.dims(2), a.dims(3), f64);
//  //array c=constant(0.0, 3, ain.dims(0), ain.dims(1), ain.dims(2), f64);
//  //std::cout << "DIMS: "<< a.dims() << "bdims: " << b.dims() << "c.dims() "
//  << c.dims() << std::endl; c(0, span, span, span)=a(1, span, span, span)*b(2,
//  span, span, span)-a(2, span, span, span)*b(1, span, span, span); c(1, span,
//  span, span)=a(2, span, span, span)*b(0, span, span, span)-a(0, span, span,
//  span)*b(2, span, span, span); c(2, span, span, span)=a(0, span, span,
//  span)*b(1, span, span, span)-a(1, span, span, span)*b(0, span, span, span);
//
//  a=reorder(a, 1, 2, 3, 0);
//  b=reorder(b, 1, 2, 3, 0);
//  c=reorder(c, 1, 2, 3, 0);
//  return c;
//};

// As well Slower
// void LLG::cross4b(array& a, array& b){ // writes result into b
//  b(span, span, span, 0)=a(span, span, span, 1)*b(span, span, span, 2)-a(span,
//  span, span, 2)*b(span, span, span, 1); b(span, span, span, 1)=a(span, span,
//  span, 2)*b(span, span, span, 0)-a(span, span, span, 0)*b(span, span, span,
//  2); b(span, span, span, 2)=a(span, span, span, 0)*b(span, span, span,
//  1)-a(span, span, span, 1)*b(span, span, span, 0);
//};
// With this in step
//    cross4b(m, heff);
//    cross4b(m, crosstemp);
//    dmdt = - state0.material.gamma/(1.+pow(state0.material.alpha, 2)) * heff -
//    state0.material.alpha*state0.material.gamma/(1.+pow(state0.material.alpha,
//    2)) * crosstemp;

// effective field H_in_Apm as sum of all terms
// std::cout << "LLG CKECK 1"<< std::endl;
// std::cout << "CHECK Exch"<< std::endl;
// array test = Exch.solve(m);
// af_print(test);
// std::cout << "CHECK Demag"<< std::endl;
// test = Demag.solve(m);
// af_print(test);
// std::cout << "xCHECK 2" << std::endl;
// std::cout << "LLG m.dims "<<m.dims()<<"heff.dims "<< heff.dims()<< std::endl;
// std::cout << "LLG CKECK 2"<< std::endl;
// std::cout << "LLG m.dims "<<m.dims()<<"heff.dims "<< heff.dims()<< std::endl;
// std::cout << "LLG CKECK 3"<< std::endl;

//// Dormand-Prince 7/8 method with stepsize control
// array LLG::DP78(const array& m, const double dt , double& err)
//{
//
//  double ptr[14], a[14][14], b[14], bhat[14], c[14];
//  ptr[1] = 0;
//  ptr[2] = 1;
//  ptr[3] = 2;
//  ptr[4] = 1;
//  ptr[5] = 3;
//  ptr[6] = 2;
//  ptr[7] = 4;
//  ptr[8] = 5;
//  ptr[9] = 6;
//  ptr[10] = 7;
//  ptr[11] = 8;
//  ptr[12] = 9;
//  ptr[13] = 1;
////C
//  a[2][1] = 5.55555555555555555555555555556e-2;
//  a[3][1] = 2.08333333333333333333333333333e-2;
//  a[3][2] = 6.25e-2;
//  a[4][1] = 3.125e-2;
//  a[4][2] = 0.e0;
//  a[4][3] = 9.375e-2;
//  a[5][1] = 3.125e-1;
//  a[5][2] = 0.e0;
//  a[5][3] = -1.171875e0;
//  a[5][4] = 1.171875e0;
//  a[6][1] = 3.75e-2;
//  a[6][2] = 0.e0;
//  a[6][3] = 0.e0;
//  a[6][4] = 1.875e-1;
//  a[6][5] = 1.5e-1;
//  a[7][1] = 4.79101371111111111111111111111e-2;
//  a[7][2] = 0.e0;
//  a[7][3] = 0.0e0;
//  a[7][4] = 1.12248712777777777777777777778e-1;
//  a[7][5] = -2.55056737777777777777777777778e-2;
//  a[7][6] = 1.28468238888888888888888888889e-2;
//  a[8][1] = 1.6917989787292281181431107136e-2;
//  a[8][2] = 0.e0;
//  a[8][3] = 0.e0;
//  a[8][4] = 3.87848278486043169526545744159e-1;
//  a[8][5] = 3.59773698515003278967008896348e-2;
//  a[8][6] = 1.96970214215666060156715256072e-1;
//  a[8][7] = -1.72713852340501838761392997002e-1;
//  a[9][1] = 6.90957533591923006485645489846e-2;
//  a[9][2] = 0.e0;
//  a[9][3] = 0.e0;
//  a[9][4] = -6.34247976728854151882807874972e-1;
//  a[9][5] = -1.61197575224604080366876923982e-1;
//  a[9][6] = 1.38650309458825255419866950133e-1;
//  a[9][7] = 9.4092861403575626972423968413e-1;
//  a[9][8] = 2.11636326481943981855372117132e-1;
//  a[10][1] = 1.83556996839045385489806023537e-1;
//  a[10][2] = 0.e0;
//  a[10][3] = 0.e0;
//  a[10][4] = -2.46876808431559245274431575997e0;
//  a[10][5] = -2.91286887816300456388002572804e-1;
//  a[10][6] = -2.6473020233117375688439799466e-2;
//  a[10][7] = 2.84783876419280044916451825422e0;
//  a[10][8] = 2.81387331469849792539403641827e-1;
//  a[10][9] = 1.23744899863314657627030212664e-1;
//  a[11][1] = -1.21542481739588805916051052503e0;
//  a[11][2] = 0.e0;
//  a[11][3] = 0.e0;
//  a[11][4] = 1.66726086659457724322804132886e1;
//  a[11][5] = 9.15741828416817960595718650451e-1;
//  a[11][6] = -6.05660580435747094755450554309e0;
//  a[11][7] = -1.60035735941561781118417064101e1;
//  a[11][8] = 1.4849303086297662557545391898e1;
//  a[11][9] = -1.33715757352898493182930413962e1;
//  a[11][10] = 5.13418264817963793317325361166e0;
//  a[12][1] = 2.58860916438264283815730932232e-1;
//  a[12][2] = 0.e0;
//  a[12][3] = 0.e0;
//  a[12][4] = -4.77448578548920511231011750971e0;
//  a[12][5] = -4.3509301377703250944070041181e-1;
//  a[12][6] = -3.04948333207224150956051286631e0;
//  a[12][7] = 5.57792003993609911742367663447e0;
//  a[12][8] = 6.15583158986104009733868912669e0;
//  a[12][9] = -5.06210458673693837007740643391e0;
//  a[12][10] = 2.19392617318067906127491429047e0;
//  a[12][11] = 1.34627998659334941535726237887e-1;
//  a[13][1] = 8.22427599626507477963168204773e-1;
//  a[13][2] = 0.e0;
//  a[13][3] = 0.e0;
//  a[13][4] = -1.16586732572776642839765530355e1;
//  a[13][5] = -7.57622116690936195881116154088e-1;
//  a[13][6] = 7.13973588159581527978269282765e-1;
//  a[13][7] = 1.20757749868900567395661704486e1;
//  a[13][8] = -2.12765911392040265639082085897e0;
//  a[13][9] = 1.99016620704895541832807169835e0;
//  a[13][10] = -2.34286471544040292660294691857e-1;
//  a[13][11] = 1.7589857770794226507310510589e-1;
//  a[13][12] = 0.e0;
////C
////C  The coefficients BHAT(*) refer to the formula used to advance the
////C  integration, here the one of order 8.  The coefficients B(*) refer
////C  to the other formula, here the one of order 7.
////C
//  bhat[1] = 4.17474911415302462220859284685e-2;
//  bhat[2] = 0.e0;
//  bhat[3] = 0.e0;
//  bhat[4] = 0.e0;
//  bhat[5] = 0.e0;
//  bhat[6] = -5.54523286112393089615218946547e-2;
//  bhat[7] = 2.39312807201180097046747354249e-1;
//  bhat[8] = 7.0351066940344302305804641089e-1;
//  bhat[9] = -7.59759613814460929884487677085e-1;
//  bhat[10] = 6.60563030922286341461378594838e-1;
//  bhat[11] = 1.58187482510123335529614838601e-1;
//  bhat[12] = -2.38109538752862804471863555306e-1;
//  bhat[13] = 2.5e-1;
////C
//  b[1] = 2.9553213676353496981964883112e-2;
//  b[2] = 0.e0;
//  b[3] = 0.e0;
//  b[4] = 0.e0;
//  b[5] = 0.e0;
//  b[6] = -8.28606276487797039766805612689e-1;
//  b[7] = 3.11240900051118327929913751627e-1;
//  b[8] = 2.46734519059988698196468570407e0;
//  b[9] = -2.54694165184190873912738007542e0;
//  b[10] = 1.44354858367677524030187495069e0;
//  b[11] = 7.94155958811272872713019541622e-2;
//  b[12] = 4.44444444444444444444444444445e-2;
//  b[13] = 0.e0;
////C
//  c[1] = 0.e0;
//  c[2] = 5.55555555555555555555555555556e-2;
//  c[3] = 8.33333333333333333333333333334e-2;
//  c[4] = 1.25e-1;
//  c[5] = 3.125e-1;
//  c[6] = 3.75e-1;
//  c[7] = 1.475e-1;
//  c[8] = 4.65e-1;
//  c[9] = 5.64865451382259575398358501426e-1;
//  c[10] = 6.5e-1;
//  c[11] = 9.24656277640504446745013574318e-1;
//  c[12] = 1.e0;
//  c[13] = c[12];
//
//  heff =  fheff(m);
//  k1   =  fdmdt(m, heff);
//  heff =  fheff(m + dt *( a[2 ][1] * k1) ); k2   =  fdmdt(m + dt *( a[2 ][1] *
//  k1) , heff); heff =  fheff(m + dt *( a[3 ][1] * k1 + a[3 ][2] * k2) ); k3 =
//  fdmdt(m + dt *( a[3 ][1] * k1 + a[3 ][2] * k2) , heff); heff =  fheff(m + dt
//  *( a[4 ][1] * k1 + a[4 ][2] * k2 +  a[4 ][3] * k3) ); k4   =  fdmdt(m + dt
//  *( a[4 ][1] * k1 + a[4 ][2] * k2 +  a[4 ][3] * k3) , heff); heff =  fheff(m
//  + dt *( a[5 ][1] * k1 + a[5 ][2] * k2 +  a[5 ][3] * k3 + a[5 ][4] * k4) );
//  k5   =  fdmdt(m + dt *( a[5 ][1] * k1 + a[5 ][2] * k2 +  a[5 ][3] * k3 + a[5
//  ][4] * k4) , heff); heff =  fheff(m + dt *( a[6 ][1] * k1 + a[6 ][2] * k2 +
//  a[6 ][3] * k3 + a[6 ][4] * k4 + a[6 ][5] * k5) ); k6   =  fdmdt(m + dt *(
//  a[6 ][1] * k1 + a[6 ][2] * k2 +  a[6 ][3] * k3 + a[6 ][4] * k4 + a[6 ][5] *
//  k5)                                                                    ,
//  heff); heff =  fheff(m + dt *( a[7 ][1] * k1 + a[7 ][2] * k2 +  a[7 ][3] *
//  k3 + a[7 ][4] * k4 + a[7 ][5] * k5 + a[7 ][6] * k6) ); k7   =  fdmdt(m + dt
//  *( a[7 ][1] * k1 + a[7 ][2] * k2 +  a[7 ][3] * k3 + a[7 ][4] * k4 + a[7 ][5]
//  * k5 + a[7 ][6] * k6)                                                    ,
//  heff); heff =  fheff(m + dt *( a[8 ][1] * k1 + a[8 ][2] * k2 +  a[8 ][3] *
//  k3 + a[8 ][4] * k4 + a[8 ][5] * k5 + a[8 ][6] * k6 + a[8 ][7] * k7) ); k8 =
//  fdmdt(m + dt *( a[8 ][1] * k1 + a[8 ][2] * k2 +  a[8 ][3] * k3 + a[8 ][4] *
//  k4 + a[8 ][5] * k5 + a[8 ][6] * k6 + a[8 ][7] * k7) , heff); heff =  fheff(m
//  + dt *( a[9 ][1] * k1 + a[9 ][2] * k2 +  a[9 ][3] * k3 + a[9 ][4] * k4 + a[9
//  ][5] * k5 + a[9 ][6] * k6 + a[9 ][7] * k7 + a[9 ][8] * k8) ); k9   = fdmdt(m
//  + dt *( a[9 ][1] * k1 + a[9 ][2] * k2 +  a[9 ][3] * k3 + a[9 ][4] * k4 + a[9
//  ][5] * k5 + a[9 ][6] * k6 + a[9 ][7] * k7 + a[9 ][8] * k8) , heff); heff =
//  fheff(m + dt *( a[10][1] * k1 + a[10][2] * k2 +  a[10][3] * k3 + a[10][4] *
//  k4 + a[10][5] * k5 + a[10][6] * k6 + a[10][7] * k7 + a[10][8] * k8 +
//  a[10][9] * k9)                              ); k10  =  fdmdt(m + dt *(
//  a[10][1] * k1 + a[10][2] * k2 +  a[10][3] * k3 + a[10][4] * k4 + a[10][5] *
//  k5 + a[10][6] * k6 + a[10][7] * k7 + a[10][8] * k8 + a[10][9] * k9) , heff);
//  heff =  fheff(m + dt *( a[11][1] * k1 + a[11][2] * k2 +  a[11][3] * k3 +
//  a[11][4] * k4 + a[11][5] * k5 + a[11][6] * k6 + a[11][7] * k7 + a[11][8] *
//  k8 + a[11][9] * k9 + a[11][10] * k10)              ); k11  =  fdmdt(m + dt
//  *( a[11][1] * k1 + a[11][2] * k2 +  a[11][3] * k3 + a[11][4] * k4 + a[11][5]
//  * k5 + a[11][6] * k6 + a[11][7] * k7 + a[11][8] * k8 + a[11][9] * k9 +
//  a[11][10] * k10)        , heff); heff =  fheff(m + dt *( a[12][1] * k1 +
//  a[12][2] * k2 +  a[12][3] * k3 + a[12][4] * k4 + a[12][5] * k5 + a[12][6] *
//  k6 + a[12][7] * k7 + a[12][8] * k8 + a[12][9] * k9 + a[12][10] * k10 +
//  a[12][11] * k11)); k12  =  fdmdt(m + dt *( a[12][1] * k1 + a[12][2] * k2 +
//  a[12][3] * k3 + a[12][4] * k4 + a[12][5] * k5 + a[12][6] * k6 + a[12][7] *
//  k7 + a[12][8] * k8 + a[12][9] * k9 + a[12][10] * k10 + a[12][11] * k11),
//  heff); heff =  fheff(m + dt *( a[13][1] * k1 + a[13][2] * k2 +  a[13][3] *
//  k3 + a[13][4] * k4 + a[13][5] * k5 + a[13][6] * k6 + a[13][7] * k7 +
//  a[13][8] * k8 + a[13][9] * k9 + a[13][10] * k10 + a[13][11] * k11 +
//  a[13][12] * k12)); k13  =  fdmdt(m + dt *( a[13][1] * k1 + a[13][2] * k2 +
//  a[13][3] * k3 + a[13][4] * k4 + a[13][5] * k5 + a[13][6] * k6 + a[13][7] *
//  k7 + a[13][8] * k8 + a[13][9] * k9 + a[13][10] * k10 + a[13][11] * k11 +
//  a[13][12] * k12), heff);
//
//  sumbk    =         dt*(bhat[1]*k1 + bhat[2]*k2+ bhat[3]*k3 + bhat[4]*k4 +
//  bhat[5]*k5 + bhat[6]*k6 + bhat[7]*k7+ bhat[8]*k8+ bhat[9]*k9 + bhat[10]*k10
//  + bhat[11]*k11 + bhat[12]*k12 + bhat[13]*k13); rk_error = sumbk -
//  dt*(b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4 + b[5]*k5 + b[6]*k6 + b[7]*k7 +
//  b[8]*k8+ b[9]*k9 + b[10]*k10 + b[11]*k11 + b[12]*k12 + b[13]*k13);
//  //rk_error = sumbk - dt*(e[1]*k1 + e[2]*k2 + e[3]*k3 + e[4]*k4 + e[5]*k5 +
//  e[6]*k6 + e[7]*k7); err=maxnorm(rk_error/givescale(max(m, m+sumbk))); return
//  sumbk;
//}

//// Bogacki 4, 5 method with stepsize control
// array LLG::BS45(const array& m, const double dt , double& err)
//{
//
//  double a[9][9];
//  double b[9];
//  double c[9];
//  double e[9];
//  a[2][1] = 1.0e0/6.0e0;
//  a[3][1] = 2.e0/27.e0;
//  a[3][2] = 4.e0/27.e0;
//  a[4][1] = 183.e0/1372.e0;
//  a[4][2] = -162.e0/343.e0;
//  a[4][3] = 1053.e0/1372.e0;
//  a[5][1] = 68.e0/297.e0;
//  a[5][2] = -4.e0/11.e0;
//  a[5][3] = 42.e0/143.e0;
//  a[5][4] = 1960.e0/3861.e0;
//  a[6][1] = 597.e0/22528.e0;
//  a[6][2] = 81.e0/352.e0;
//  a[6][3] = 63099.e0/585728.e0;
//  a[6][4] = 58653.e0/366080.e0;
//  a[6][5] = 4617.e0/20480.e0;
//  a[7][1] = 174197.e0/959244.e0;
//  a[7][2] = -30942.e0/79937.e0;
//  a[7][3] = 8152137.e0/19744439.e0;
//  a[7][4] = 666106.e0/1039181.e0;
//  a[7][5] = -29421.e0/29068.e0;
//  a[7][6] = 482048.e0/414219.e0;
//  a[8][1] = 587.e0/8064.e0;
//  a[8][2] = 0.e0;
//  a[8][3] = 4440339.e0/15491840.e0;
//  a[8][4] = 24353.e0/124800.e0;
//  a[8][5] = 387.e0/44800.e0;
//  a[8][6] = 2152.e0/5985.e0;
//  a[8][7] = 7267.e0/94080.e0;
////C  The coefficients B(*) refer to the formula of order 4.
////C
//  b[1] = 2479.e0/34992.e0;
//  b[2] = 0.e0;
//  b[3] = 123.e0/416.e0;
//  b[4] = 612941.e0/3411720.e0;
//  b[5] = 43.e0/1440.e0;
//  b[6] = 2272.e0/6561.e0;
//  b[7] = 79937.e0/1113912.e0;
//  b[8] = 3293.e0/556956.e0;
////C  The coefficients E(*) refer to an estimate of the local error based on
////C  the first formula of order 4.  It is the difference of the fifth order
////C  result, here located in A(8, *), and the fourth order result.  By
////C  construction both E(2) and E(7) are zero.
////C
//  e[1] = -3.e0/1280.e0;
//  e[2] = 0.e0;
//  e[3] = 6561.e0/632320.e0;
//  e[4] = -343.e0/20800.e0;
//  e[5] = 243.e0/12800.e0;
//  e[6] = -1.e0/95.e0;
//  e[7] = 0.e0;
////C
//  c[1] = 0.e0;
//  c[2] = 1.e0/6.e0;
//  c[3] = 2.e0/9.e0;
//  c[4] = 3.e0/7.e0;
//  c[5] = 2.e0/3.e0;
//  c[6] = 3.e0/4.e0;
//  c[7] = 1.e0;
//  c[8] = 1.e0;
//
//  int s=8;
//  if(reject || calls==0 || llg_wasnormalized){
//    heff =  fheff(m);
//    k[1]   =  fdmdt(m, heff);
//    }
//  else
//    k[1]=k[s];
//  for(int i=2;i<=s;i++){
//      rktemp=constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3,
//      f64);
//    for(int j=1;j<i;j++){
//      //std::cout<<"a"<<i<<j<<std::endl;
//      rktemp+=a[i][j] * k[j];
//    }
//    rktemp*=dt;
//    heff= fheff(m + rktemp);
//    k[i]= fdmdt(m + rktemp, heff);
//  //std::cout<<std::endl;
//  }
//  sumbk=constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3, f64);
//  for(int i=1;i<s;i++){
//    sumbk+=a[s][i]*k[i];
//  }
//  sumbk*=dt;
//
//  //rk_error=constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3,
//  f64);
//  //for(int i=1;i<=s;i++){
//  //  rk_error+=b[i]*k[i];
//  //}
//  //rk_error*=dt;
//  //rk_error=sumbk-rk_error;
//
//  rk_error=constant(0.0, state0.mesh.nx, state0.mesh.ny, state0.mesh.nz, 3,
//  f64); for(int i=1;i<s;i++){//i=1, 2, ..., 7
//    rk_error+=e[i]*k[i];
//  }
//  rk_error*=dt;
//  err=maxnorm(rk_error/givescale(max(m, m+sumbk)));
//  //rk_error=sumbk-rk_error;
//
//  //TODO rk_error = sumbk - dt*(e[1]*k1 + e[2]*k2 + e[3]*k3 + e[4]*k4 +
//  e[5]*k5 + e[6]*k6 + e[7]*k7 + ?);
//
////  std::cout<<"h of k1 = "<<afvalue((heff)(0, 0, 0, 0))<<"\n"<<std::endl;
////  std::cout<<"k1 = "<<afvalue((k1)(0, 0, 0, 0))<<"\n"<<std::endl;
////  std::cout<<"m k1 = "<<afvalue((m)(0, 0, 0, 0))<<"\t"<<"h of k1 =
///"<<afvalue((heff)(0, 0, 0, 0))<<"\t"<<"k1= "<<afvalue(k1(0, 0, 0, 0))<<"\t
/// dtfdmdt= "<<afvalue(fdmdt(m, heff)(0, 0, 0, 0))<<"\n"<<std::endl; / Start
/// Alternative /  if(reject || calls==0 || llg_wasnormalized){ /    heff =
/// fheff(m); /    k1   =  fdmdt(m, heff); /    } /  else /    k1=k8; /  heff =
/// fheff(m + dt *( a[2][1] * k1) ); /  k2   =  fdmdt(m + dt *( a[2][1] * k1) ,
/// heff); /  heff =  fheff(m + dt *( a[3][1] * k1 + a[3][2] * k2) ); /  k3   =
/// fdmdt(m + dt *( a[3][1] * k1 + a[3][2] * k2) , heff); /  heff =  fheff(m +
/// dt
///*( a[4][1] * k1 + a[4][2] * k2 +  a[4][3] * k3) ); /  k4   =  fdmdt(m + dt *(
/// a[4][1] * k1 + a[4][2] * k2 +  a[4][3] * k3) , heff); /  heff =  fheff(m +
/// dt
///*( a[5][1] * k1 + a[5][2] * k2 +  a[5][3] * k3 + a[5][4] * k4) ); /  k5   =
/// fdmdt(m + dt *( a[5][1] * k1 + a[5][2] * k2 +  a[5][3] * k3 + a[5][4] * k4)
///, heff); /  heff =  fheff(m + dt *( a[6][1] * k1 + a[6][2] * k2 +  a[6][3] *
/// k3 + a[6][4] * k4 + a[6][5] * k5)                                    ); / k6
///=  fdmdt(m + dt *( a[6][1] * k1 + a[6][2] * k2 +  a[6][3] * k3 + a[6][4] * k4
///+ a[6][5] * k5)                              , heff); /  heff =  fheff(m + dt
///*( a[7][1] * k1 + a[7][2] * k2 +  a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5
///+ a[7][6] * k6)                     ); /  k7   =  fdmdt(m + dt *( a[7][1] *
/// k1 + a[7][2] * k2 +  a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] *
/// k6)               , heff); /  heff =  fheff(m + dt *( a[8][1] * k1 + a[8][2]
///* k2 +  a[8][3] * k3 + a[8][4] * k4 + a[8][5] * k5 + a[8][6] * k6 + a[8][7] *
/// k7)      ); /  k8   =  fdmdt(m + dt *( a[8][1] * k1 + a[8][2] * k2 + a[8][3]
///* k3 + a[8][4] * k4 + a[8][5] * k5 + a[8][6] * k6 + a[8][7] * k7), heff);
////
////  sumbk    =         dt*(a[8][1]*k1 + a[8][2]*k2+ a[8][3]*k3 + a[8][4]*k4 +
/// a[8][5]*k5 + a[8][6]*k6 + a[8][7]*k7); /  rk_error = sumbk - dt*(b[1]*k1 +
/// b[2]*k2 + b[3]*k3 + b[4]*k4 + b[5]*k5 + b[6]*k6 + b[7]*k7 + b[8]*k8);
//  //TODO rk_error = sumbk - dt*(e[1]*k1 + e[2]*k2 + e[3]*k3 + e[4]*k4 +
//  e[5]*k5 + e[6]*k6 + e[7]*k7);
//// End Alternative
//  return sumbk;
//}
//// Tsitorous 4/5
// array LLG::tsit45(const array& m, const double dt, double& err)
//{
//  //state0.t+=dt;
//  const double c2 = 0.161;
//  const double c3 = 0.327;
//  const double c4 = 0.9;
//  const double c5 = 0.9800255409045097;
//  const double a32 = 0.3354806554923570;
//  const double a42 = -6.359448489975075;
//  const double a43 = 4.362295432869581;
//  const double a52 = -11.74888356406283;
//  const double a53 = 7.495539342889836;
//  const double a54 = -0.09249506636175525;
//  const double a62 = -12.92096931784711;
//  const double a64 = -0.07158497328140100;
//  const double a63 = 8.159367898576159;
//  const double a65 = -0.02826905039406838;
//  const double b1 = 0.09646076681806523;
//  const double b3 = 0.4798896504144996;
//  const double b5 = -3.290069515436081;
//  const double b2 = 0.01;
//  const double b4 = 1.379008574103742;
//  const double b6 = 2.324710524099774;
//  const double b1_hat = 0.001780011052226;
//  const double b2_hat = 0.000816434459657;
//  const double b3_hat = -0.007880878010262;
//  const double b4_hat = 0.144711007173263;
//  const double b5_hat = -0.582357165452555;
//  const double b6_hat = 0.458082105929187;
//  const double b7_hat = 1./66.;
//  //double a7i = b i , i = 1, 2, · · · , 6
//  //
//  //double c[7]={c1, c2, c3, c4, c5, c6, c7};
//  const double b[8]={0, b1, b2, b3, b4, b5, b6, 0};
//  const double b_hat[8]={0, b1_hat, b2_hat, b3_hat, b4_hat, b5_hat, b6_hat,
//  b7_hat}; const double a[8][8]={
//    {0, 0, 0, 0, 0, 0, 0, 0},
//    {0, 0, 0, 0, 0, 0, 0, 0},
//    {0, c2, 0, 0, 0, 0, 0, 0},
//    {0, c3-a32, a32, 0, 0, 0, 0, 0},
//    {0, c4 - a42 - a43, a42, a43, 0, 0, 0, 0},
//    {0, c5 - a52 - a53 - a54, a52, a53, a54, 0, 0, 0},
//    {0, 1 - a62 - a63 - a64 - a65, a62, a63, a64, a65, 0, 0},
//    {0, b1, b2, b3, b4, b5, b6, 0.}
//  };
//
//  heff =       fheff(m ); k1   =  dt * fdmdt(m , heff); heff =       fheff(m
//  + a[2][1] * k1 ); k2   =  dt * fdmdt(m   + a[2][1] * k1 , heff); heff =
//  fheff(m   + a[3][1] * k1 + a[3][2]  * k2 ); k3   =  dt * fdmdt(m   + a[3][1]
//  * k1 + a[3][2]  * k2 , heff); heff =       fheff(m   + a[4][1] * k1 +
//  a[4][2]  * k2 +  a[4][3]  * k3 ); k4   =  dt * fdmdt(m   + a[4][1] * k1 +
//  a[4][2]  * k2 +  a[4][3]  * k3                                             ,
//  heff); heff =       fheff(m   + a[5][1] * k1 + a[5][2]  * k2 +  a[5][3]  *
//  k3 + a[5][4] * k4                                    ); k5   =  dt * fdmdt(m
//  + a[5][1] * k1 + a[5][2]  * k2 +  a[5][3]  * k3 + a[5][4] * k4 , heff); heff
//  =       fheff(m   + a[6][1] * k1 + a[6][2]  * k2 +  a[6][3]  * k3 + a[6][4]
//  * k4 + a[6][5] * k5                     ); k6   =  dt * fdmdt(m   + a[6][1]
//  * k1 + a[6][2]  * k2 +  a[6][3]  * k3 + a[6][4] * k4 + a[6][5] * k5 , heff);
//  heff =       fheff(m   + a[7][1] * k1 + a[7][2]  * k2 +  a[7][3]  * k3 +
//  a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6      ); k7   =  dt * fdmdt(m   +
//  a[7][1] * k1 + a[7][2]  * k2 +  a[7][3]  * k3 + a[7][4] * k4 + a[7][5] * k5
//  + a[7][6] * k6, heff);
//
//  sumbk    =  b[1] * k1 + b[2]  * k2 +  b[3]  * k3 + b[4] * k4 + b[5] * k5 +
//  b[6] * k6;
//  //rk_error = (  b_hat[1] * k1 + b_hat[2]  * k2 +  b_hat[3]  * k3 + b_hat[4]
//  * k4 + b_hat[5] * k5 + b_hat[6] * k6 + b_hat[7] * k7 ); rk_error = sumbk - (
//  b_hat[1] * k1 + b_hat[2]  * k2 +  b_hat[3]  * k3 + b_hat[4] * k4 + b_hat[5]
//  * k5 + b_hat[6] * k6 + b_hat[7] * k7 );
//
//  err=maxnorm(rk_error/givescale(max(m, m+sumbk)));
//  //err = maxnorm(rk_error);
//  return sumbk;
//}

//// Dormand-Prince 4/5 method
// array LLG::DP45(const array& m, const double dt, double& err)
//{
//  //double c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;
//  double a21=0.2, a31=3.0/40.0, a32=9.0/40.0, a41=44.0/45.0, a42=-56.0/15.0,
//  a43=32.0/9.0, a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0,
//  a54=-212.0/729.0, a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0,
//  a64=49.0/176.0, a65=-5103.0/18656.0, a71=35.0/384.0, a73=500.0/1113.0,
//  a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0, e1=71.0/57600.0,
//  e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0, e6=22.0/525.0,
//  e7=-1.0/40.0; const double e[8]={0, e1, 0, e3, e4, e5, e6, e7}; const double
//  a[8][8]={
//    {0, 0, 0, 0, 0, 0, 0, 0},
//    {0, 0, 0, 0, 0, 0, 0, 0},
//    {0, a21, 0, 0, 0, 0, 0, 0},
//    {0, a31, a32, 0, 0, 0, 0, 0},
//    {0, a41, a42, a43, 0, 0, 0, 0},
//    {0, a51, a52, a53, a54, 0, 0, 0},
//    {0, a61, a62, a63, a64, a65, 0, 0},
//    {0, a71, 0  , a73, a74, a75, a76, 0.}
//  };
//
////  std::cout<<"m of k1 = "<<afvalue((m)(0, 0, 0, 0))<<"\n"<<std::endl;
////  std::cout<<"h of k1 = "<<afvalue((heff)(0, 0, 0, 0))<<std::endl;
////  std::cout<<"k1= "<<afvalue(k1(0, 0, 0, 0))<<"\t"<<afvalue(k1(1, 1, 0,
/// 1))<<"\n"<<std::endl;
//  if(reject) std::cout<<"!!!!!!!! Prev was rejected"<<std::endl; if(reject){
//    heff =       fheff(m                               );
//    k1   =  dt * fdmdt(m                         , heff);
//  }
//  else
//    k1=k7;
//
//  mtemp= m + a[2][1] * k1;
//  heff =       fheff(mtemp                           );
//  k2   =  dt * fdmdt(mtemp                     , heff);
//
//  mtemp=m + a[3][1] * k1 + a[3][2] * k2;
//  heff =       fheff(mtemp                           );
//  k3   =  dt * fdmdt(mtemp                     , heff);
//
//  mtemp=m + a[4][1] * k1 + a[4][2] * k2 +  a[4][3] * k3;
//  heff =       fheff(mtemp                           );
//  k4   =  dt * fdmdt(mtemp                     , heff);
//
//  mtemp=m + a[5][1] * k1 + a[5][2] * k2 +  a[5][3] * k3 + a[5][4] * k4;
//  heff =       fheff(mtemp                                    );
//  k5   =  dt * fdmdt(mtemp                              , heff);
//
//  mtemp=m + a[6][1] * k1 + a[6][2] * k2 +  a[6][3] * k3 + a[6][4] * k4 +
//  a[6][5] * k5; heff =       fheff(mtemp                     ); k6   =  dt *
//  fdmdt(mtemp               , heff);
//
//  sumbk=    a[7][1] * k1                 + a[7][3] * k3 + a[7][4] * k4 +
//  a[7][5] * k5 + a[7][6] * k6; heff =       fheff(m + sumbk      ); k7   =  dt
//  * fdmdt(m + sumbk, heff);
////  std::cout<<"m of k7 = "<<afvalue((m + a[7][1]*k1 + a[7][3]*k3 + a[7][4]*k4
///+ a[7][5]*k5 + a[7][6]*k6)(0, 0, 0, 0))<<std::endl; /  std::cout<<"h of k7 =
///"<<afvalue((heff)(0, 0, 0, 0))<<std::endl; /  std::cout<<"k7=
///"<<afvalue(k7(0, 0, 0, 0))<<"\t"<<afvalue(k7(1, 1, 0, 1))<<std::endl;
//
//  //Todo: not differs in 7th digit
//  //std::cout.precision(12);
//  //std::cout << maxnorm(k1(0, 0, 0, 0)) << "\t" << maxnorm(k7(0, 0, 0, 0)) <<
//  std::endl;
//
////  std::cout<<"m+sumbk = "<<afvalue((m+sumbk)(0, 0, 0, 0))<<std::endl;
//  //rk_error = (a[7][1] - e[1] ) * k1 + (a[7][2] - e[2] * k2 ) + (a[7][3] -
//  e[3] ) * k3 + ( a[7][4] -e[4] ) * k4 + (a[7][5] - e[5]) * k5 + (a[7][6] -
//  e[6] )* k6 -e[7] * k7;
//
//  rk_error = sumbk - (e[1]*k1 + e[2]*k2 + e[3]*k3 + e[4]*k4 + e[5]*k5 +
//  e[6]*k6 + e[7]*k7); err=maxnorm(rk_error/givescale(max(m, m+sumbk)));
//  //std::cout << "h = "<<h<< " error = "<<err<<" maxnorm(rk_error)
//  "<<maxnorm(rk_error)<<" givescale "<< maxnorm(givescale(max(m,
//  m+sumbk)))<<std::endl;
//  //rk_abs_error = maxnorm(rk_error);
//  return sumbk;
//}
