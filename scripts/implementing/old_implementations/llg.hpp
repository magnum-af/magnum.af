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
    float E(const State& state);
    float cpu_time();
    void print_cpu_time(std::ostream&);
    // vtk fieldterm writer
    void write_fieldterms_atom (const State& state, const std::string filepath);
    void write_fieldterms_micro(const State& state, const std::string filepath);
    void relax(State& state, float precision = 1e-10, const int iloop = 100, const int iwritecout = 1000);

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
    float llg_normtol{1e-5};// If maxnorm(vecnorm(mtemp)) exceedes this value, renormalize is called
    int  llg_normalize_counter{0};// number of normalizations performed

    // Counting step calls
    int calls{0};

    // Computation Arrays
    array heff, crosstemp, dmdt, sumbk, rk_error, mtemp, rktemp;

    // Explicit Euler Integrator
    array explicitEuler(const array&, float dt);

    // RK-methods
    array k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
    array k[14];
    int s=0;
    bool FSAL=false;
    float a[14][14]={{0.}};
    float b[14]={0.};
    float bhat[14]={0.};
    float c[14]={0.};
    float e[14]={0.};

    array rk4(const array& , const float dt);
    array rk4minimal(const array&, const float dt);
    array rk4_3o8(const array&, const float dt); //RK4-3/8 method


    //RKF
    //float rk_abs_error{0.};
    //float rk_abs_tol_error{8.e-6}; // 8.e-6 is pretty good
    //float rk_rel_error{0.};
    //float rk_rel_tol_error{5.e-3}; // 5.e-3 is pretty good
    //float h_abs{0.};
    //float h_rel{0.};

    float h; // working variable in Controller, actual stepsize and next stepsize, min of h_abs and h_rel
    float h_stepped; // step size set after step is called
    //int counter_abs{0};
    //int counter_rel{0};
    //int counter_abs_reject{0};
    //int counter_rel_reject{0};
    //int p{5};
    // Runge-Kutta-Fehlberg
    array RKF5(const array&, const float dt);

    // Runge-Kutta-Fehlberg with stepsize control
    array RKF45(const array&, const float dt, float& err);

    // Cash-Karp  with stepsize control
    array CK45(const array&, const float dt, float& err);

    // Tsitorous 4/5  with stepsize control
    array tsit45(const array&, const float dt, float& err);

    // Dormand-Prince  with stepsize control
    array DP45(const array&, const float dt, float& err);
    array DP78(const array&, const float dt, float& err);
    //array DP45(array m, float dt, float& rk_abs_error);// Dormand-Prince  with stepsize control

    // Bogacki-Shampine with stepsize control
    array BS45(const array&, const float dt, float& err);
    array BS45de(const array&, const float dt, float& err);//Double error estimation
    array BS23(const array&, const float h, float& err);

    //Numerical Recipies Adaptive Stepsize control
    //const float err0{1.};// Desired error
    //const float atol{1e-8};//Tolerated absolute error
    //const float rtol{1e-8};//Tolerated relative error
    //const float hmin{1e-15};
    //const float hmax{3.5e-10};
    array  givescale(const array& a); // Scale function return= atol + abs(y) * rtol
    float  err{.0};      // Estimated error
    int counter_reject{0};// # of rejections
    int counter_accepted{0};// # of accepced steps
    int counter_hmax{0};// # of rejections
    int counter_hmin{0};// # of rejections
    //For Controller
    //bool controller(const float err, float& h);
    Controller controller = Controller();
    bool reject{false};
    float errold{1.0e-4};
    int counter_maxscale{0};// # of rejections
    int counter_minscale{0};// # of rejections
    float hnext;

    // Timing
    float time_integrator{0.};
    timer timer_integrator;
    float time_heff{0.};
    timer timer_heff;
    float time_fdmdt{0.};
    timer timer_fdmdt;
};
#endif
