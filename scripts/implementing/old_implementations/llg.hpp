#ifndef LLG_H
#define LLG_H
#include <vector>
#include <memory>
#include <initializer_list>
#include "arrayfire.h"
#include "../llg_terms/LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
#include "controller.hpp"

using namespace af;

class LLG {
  public:
    LLG (State, std::vector<std::shared_ptr<LLGTerm> >);
    //Field Term Contributions
    std::vector<std::shared_ptr<LLGTerm> > Fieldterms;
    // Class Objects
    State state0;

    //Energy calculation
    double E(const State& state);
    double cpu_time();
    void print_cpu_time(std::ostream&);
    // vtk fieldterm writer
    void write_fieldterms_atom (const State& state, const std::string filepath);
    void write_fieldterms_micro(const State& state, const std::string filepath);
    void relax(State& state, double precision = 1e-10, const int iloop = 100, const int iwritecout = 1000);

    // rhs of the LLG equation
    bool fdmdt_dissipation_term_only{false};//If ture, only use the energy dissipation term (i.e. Mx(MxH)) in the LLG
    array fdmdt(const array& m, const array& heff);
    int calls_fdmdt{0};
    array fdmdtminimal(array m, array heff);

    // Calculation of effective field with optional zeeman field
    array fheff(const array& m);
    array fheffminimal(array m);
    long int h_addr(const State& state);//Getter function for effective field for wrapping

    std::vector<array> h_addr_temp_array;
    //Alternative:
    //array h_addr_temp_array;

    array step(State& state);
    bool llg_wasnormalized{true};//set true after normalization in step, false if not normalized
                                 // DP45 uses this to decide wether k1 and heff has to be calculated
                                 // -> leave the value true if normalization is on default
    double llg_normtol{1e-5};// If maxnorm(vecnorm(mtemp)) exceedes this value, renormalize is called
    int  llg_normalize_counter{0};// number of normalizations performed

    // Counting step calls
    int calls{0};

    // Computation Arrays
    array heff, crosstemp, dmdt, sumbk, rk_error, mtemp, rktemp;

    // Explicit Euler Integrator
    array explicitEuler(const array&, double dt);

    // RK-methods
    array k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
    array k[14];
    int s=0;
    bool FSAL=false;
    double a[14][14]={{0.}};
    double b[14]={0.};
    double bhat[14]={0.};
    double c[14]={0.};
    double e[14]={0.};

    array rk4(const array& , const double dt);
    array rk4minimal(const array&, const double dt);
    array rk4_3o8(const array&, const double dt); //RK4-3/8 method


    //RKF
    //double rk_abs_error{0.};
    //double rk_abs_tol_error{8.e-6}; // 8.e-6 is pretty good
    //double rk_rel_error{0.};
    //double rk_rel_tol_error{5.e-3}; // 5.e-3 is pretty good
    //double h_abs{0.};
    //double h_rel{0.};

    double h; // working variable in Controller, actual stepsize and next stepsize, min of h_abs and h_rel
    double h_stepped; // step size set after step is called
    //int counter_abs{0};
    //int counter_rel{0};
    //int counter_abs_reject{0};
    //int counter_rel_reject{0};
    //int p{5};
    // Runge-Kutta-Fehlberg
    array RKF5(const array&, const double dt);

    // Runge-Kutta-Fehlberg with stepsize control
    array RKF45(const array&, const double dt, double& err);

    // Cash-Karp  with stepsize control
    array CK45(const array&, const double dt, double& err);

    // Tsitorous 4/5  with stepsize control
    array tsit45(const array&, const double dt, double& err);

    // Dormand-Prince  with stepsize control
    array DP45(const array&, const double dt, double& err);
    array DP78(const array&, const double dt, double& err);
    //array DP45(array m, double dt, double& rk_abs_error);// Dormand-Prince  with stepsize control

    // Bogacki-Shampine with stepsize control
    array BS45(const array&, const double dt, double& err);
    array BS45de(const array&, const double dt, double& err);//Double error estimation
    array BS23(const array&, const double h, double& err);

    //Numerical Recipies Adaptive Stepsize control
    //const double err0{1.};// Desired error
    //const double atol{1e-8};//Tolerated absolute error
    //const double rtol{1e-8};//Tolerated relative error
    //const double hmin{1e-15};
    //const double hmax{3.5e-10};
    array  givescale(const array& a); // Scale function return= atol + abs(y) * rtol
    double  err{.0};      // Estimated error
    int counter_reject{0};// # of rejections
    int counter_accepted{0};// # of accepced steps
    int counter_hmax{0};// # of rejections
    int counter_hmin{0};// # of rejections
    //For Controller
    //bool controller(const double err, double& h);
    Controller controller = Controller();
    bool reject{false};
    double errold{1.0e-4};
    int counter_maxscale{0};// # of rejections
    int counter_minscale{0};// # of rejections
    double hnext;

    // Timing
    double time_integrator{0.};
    timer timer_integrator;
    double time_heff{0.};
    timer timer_heff;
    double time_fdmdt{0.};
    timer timer_fdmdt;
};
#endif
