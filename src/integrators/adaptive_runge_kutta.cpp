#include "adaptive_runge_kutta.hpp"
#include "util/func.hpp"
#include "util/misc.hpp"

namespace magnumafcpp {

AdaptiveRungeKutta::AdaptiveRungeKutta(std::string scheme_, Controller controller_, const bool normalize_,
                                       const bool verbose)
    : scheme_(scheme_), controller_(controller_), normalize_(normalize_) {
    if (scheme_ == "RKF45") {
        if (verbose)
            printf("Adaptive Runge Kutta: Initializing RKF45 method.\n");
    } else if (scheme_ == "DP45") {
        if (verbose)
            printf("Adaptive Runge Kutta: Initializing DP45 method.\n");
    } else if (scheme_ == "BS45") {
        if (verbose)
            printf("Adaptive Runge Kutta: Initializing BS45 method.\n");
    } else if (scheme_ == "DP78") {
        if (verbose)
            printf("Adaptive Runge Kutta: Initializing DP78 method.\n");
    } else if (scheme_ == "BS23") {
        if (verbose)
            printf("Adaptive Runge Kutta: Initializing BS23 method.\n");
    } else {
        printf("%s Integration method not found. Valid arguments are 'RKF45' "
               "(default) , 'DP45', 'BS45', 'DP78', 'BS23'\n",
               bold_red("Error:").c_str());
        exit(EXIT_FAILURE);
    }
}

void AdaptiveRungeKutta::step(State& state) {
    af::timer timer_allsteps = af::timer::start();
    af::array mtemp;
    do {
        if (scheme_ == "RKF45") {
            mtemp = RKF45(state, h_, err_);
        } else if (scheme_ == "DP45") {
            mtemp = DP45(state, h_, err_);
        } else if (scheme_ == "BS45") {
            mtemp = BS45(state, h_, err_);
        } else if (scheme_ == "DP78") {
            mtemp = DP78(state, h_, err_);
        } else {
            mtemp = BS23(state, h_, err_);
        }
    } while (!controller_.success(err_, h_));

    state.t += h_; // h is the actual timestep taken by the controller_
    h_ = controller_.get_hnext();
    state.m += mtemp;
    if (normalize_) {
        if (state.Ms_field.isempty()) {
            state.m = normalize(state.m);
        } else {
            state.m = normalize_handle_zero_vectors(state.m);
        }
    }
    time_allsteps_ += af::timer::stop(timer_allsteps);
    state.steps++;
    this->accumulated_steps++;
    af::eval(state.m);
}

// Runge-Kutta-Fehlberg method with stepsize control
af::array AdaptiveRungeKutta::RKF45(const State& state, const double dt, double& err_) {
    State tempstate = state;
    // stage1
    af::array k1 = dt * f(state);

    // stage2
    tempstate.t = state.t + 1. / 4. * dt;
    tempstate.m = state.m + 1. / 4. * k1;
    af::array k2 = dt * f(tempstate);

    // stage3
    tempstate.t = state.t + 3. / 8. * dt;
    tempstate.m = state.m + 3. / 32. * k1 + 9 / 32. * k2;
    af::array k3 = dt * f(tempstate);

    // stage4
    tempstate.t = state.t + 12. / 13. * dt;
    tempstate.m = state.m + 1932. / 2197. * k1 - 7200. / 2197. * k2 + 7296. / 2197. * k3;
    af::array k4 = dt * f(tempstate);

    // stage5
    tempstate.t = state.t + dt;
    tempstate.m = state.m + 439. / 216. * k1 - 8. * k2 + 3680. / 513. * k3 - 845. / 4104. * k4;
    af::array k5 = dt * f(tempstate);

    // stage6
    tempstate.t = state.t + 1. / 2. * dt;
    tempstate.m = state.m - 8. / 27. * k1 + 2. * k2 - 3544. / 2565. * k3 + 1859. / 4104. * k4 - 11. / 40. * k5;
    af::array k6 = dt * f(tempstate);

    af::array sumbk = 16. / 135. * k1 + 6656. / 12825. * k3 + 28561. / 56430. * k4 - 9. / 50. * k5 + 2. / 55. * k6;
    af::array rk_error = sumbk - (25. / 216. * k1 + 1408. / 2565. * k3 + 2197. / 4104. * k4 - 1. / 5. * k5);

    err_ = max_4d_abs(rk_error / controller_.givescale(max(state.m, state.m + sumbk)));
    return sumbk;
}

// Dormand-Prince 4/5 method
af::array AdaptiveRungeKutta::DP45(const State& state, const double dt, double& err_) {
    State tempstate = state;

    double a[8][7] = {{0}};
    double e[8] = {0};
    double c[8] = {0};

    c[2] = 0.2, c[3] = 0.3, c[4] = 0.8, c[5] = 8.0 / 9.0, c[6] = 1, c[7] = 1;
    e[1] = 71.0 / 57600.0, e[3] = -71.0 / 16695.0, e[4] = 71.0 / 1920.0, e[5] = -17253.0 / 339200.0,
    e[6] = 22.0 / 525.0, e[7] = -1.0 / 40.0, a[2][1] = 0.2, a[3][1] = 3.0 / 40.0, a[3][2] = 9.0 / 40.0,
    a[4][1] = 44.0 / 45.0, a[4][2] = -56.0 / 15.0, a[4][3] = 32.0 / 9.0, a[5][1] = 19372.0 / 6561.0,
    a[5][2] = -25360.0 / 2187.0, a[5][3] = 64448.0 / 6561.0, a[5][4] = -212.0 / 729.0, a[6][1] = 9017.0 / 3168.0,
    a[6][2] = -355.0 / 33.0, a[6][3] = 46732.0 / 5247.0, a[6][4] = 49.0 / 176.0, a[6][5] = -5103.0 / 18656.0,
    a[7][1] = 35.0 / 384.0, a[7][3] = 500.0 / 1113.0, a[7][4] = 125.0 / 192.0, a[7][5] = -2187.0 / 6784.0,
    a[7][6] = 11.0 / 84.0;

    // Stage 1
    af::array k1;
    if (controller_.get_reject() || normalize_ || state.steps == 0) {
        k1 = dt * f(tempstate);
    } else {
        k1 = k_FSAL;
    }

    // Stage 2
    tempstate.t = state.t + c[2] * dt;
    tempstate.m = state.m + a[2][1] * k1;
    af::array k2 = dt * f(tempstate);

    // Stage 3
    tempstate.t = state.t + c[3] * dt;
    tempstate.m = state.m + a[3][1] * k1 + a[3][2] * k2;
    af::array k3 = dt * f(tempstate);

    // Stage 4
    tempstate.t = state.t + c[4] * dt;
    tempstate.m = state.m + a[4][1] * k1 + a[4][2] * k2 + a[4][3] * k3;
    af::array k4 = dt * f(tempstate);

    // Stage 5
    tempstate.t = state.t + c[5] * dt;
    tempstate.m = state.m + a[5][1] * k1 + a[5][2] * k2 + a[5][3] * k3 + a[5][4] * k4;
    af::array k5 = dt * f(tempstate);

    // Stage 6
    tempstate.t = state.t + c[6] * dt;
    tempstate.m = state.m + a[6][1] * k1 + a[6][2] * k2 + a[6][3] * k3 + a[6][4] * k4 + a[6][5] * k5;
    af::array k6 = dt * f(tempstate);

    // Stage 7
    tempstate.t = state.t + c[7] * dt;
    tempstate.m = state.m + a[7][1] * k1 + a[7][2] * k2 + a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6;
    k_FSAL = dt * f(tempstate);

    af::array sumbk = a[7][1] * k1 + a[7][2] * k2 + a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6;
    af::array rk_error = e[1] * k1 + e[2] * k2 + e[3] * k3 + e[4] * k4 + e[5] * k5 + e[6] * k6 + e[7] * k_FSAL;

    err_ = max_4d_abs(rk_error / controller_.givescale(max(state.m, state.m + sumbk)));
    return sumbk;
}

// Bogacki 4, 5 method with sigle error andstepsize control
af::array AdaptiveRungeKutta::BS45(const State& state, const double dt, double& err_) {
    State tempstate = state;

    double a[9][8] = {{0.}};
    double b[9] = {0.};
    double c[9] = {0.};

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

    c[1] = 0.e0;
    c[2] = 1.e0 / 6.e0;
    c[3] = 2.e0 / 9.e0;
    c[4] = 3.e0 / 7.e0;
    c[5] = 2.e0 / 3.e0;
    c[6] = 3.e0 / 4.e0;
    c[7] = 1.e0;
    c[8] = 1.e0;

    b[1] = 2479.e0 / 34992.e0;
    b[2] = 0.e0;
    b[3] = 123.e0 / 416.e0;
    b[4] = 612941.e0 / 3411720.e0;
    b[5] = 43.e0 / 1440.e0;
    b[6] = 2272.e0 / 6561.e0;
    b[7] = 79937.e0 / 1113912.e0;
    b[8] = 3293.e0 / 556956.e0;

    // Stage 1
    af::array k1;
    if (controller_.get_reject() || normalize_ || step_calls_ == 0) {
        k1 = dt * f(tempstate);
    } else {
        k1 = k_FSAL;
    }

    // Stage 2
    tempstate.t = state.t + c[2] * dt;
    tempstate.m = state.m + a[2][1] * k1;
    af::array k2 = dt * f(tempstate);

    // Stage 3
    tempstate.t = state.t + c[3] * dt;
    tempstate.m = state.m + a[3][1] * k1 + a[3][2] * k2;
    af::array k3 = dt * f(tempstate);

    // Stage 4
    tempstate.t = state.t + c[4] * dt;
    tempstate.m = state.m + a[4][1] * k1 + a[4][2] * k2 + a[4][3] * k3;
    af::array k4 = dt * f(tempstate);

    // Stage 5
    tempstate.t = state.t + c[5] * dt;
    tempstate.m = state.m + a[5][1] * k1 + a[5][2] * k2 + a[5][3] * k3 + a[5][4] * k4;
    af::array k5 = dt * f(tempstate);

    // Stage 6
    tempstate.t = state.t + c[6] * dt;
    tempstate.m = state.m + a[6][1] * k1 + a[6][2] * k2 + a[6][3] * k3 + a[6][4] * k4 + a[6][5] * k5;
    af::array k6 = dt * f(tempstate);

    // Stage 7
    tempstate.t = state.t + c[7] * dt;
    tempstate.m = state.m + a[7][1] * k1 + a[7][2] * k2 + a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6;
    af::array k7 = dt * f(tempstate);

    // Stage 8
    tempstate.t = state.t + c[8] * dt;
    tempstate.m = state.m + a[8][1] * k1 + a[8][2] * k2 + a[8][3] * k3 + a[8][4] * k4 + a[8][5] * k5 + a[8][6] * k6 +
                  a[8][7] * k7;
    k_FSAL = dt * f(tempstate);

    af::array sumbk =
        a[8][1] * k1 + a[8][2] * k2 + a[8][3] * k3 + a[8][4] * k4 + a[8][5] * k5 + a[8][6] * k6 + a[8][7] * k7;
    af::array rk_error =
        sumbk - (b[1] * k1 + b[2] * k2 + b[3] * k3 + b[4] * k4 + b[5] * k5 + b[6] * k6 + b[7] * k7 + b[8] * k_FSAL);
    err_ = max_4d_abs(rk_error / controller_.givescale(max(state.m, state.m + sumbk)));
    return sumbk;
}

// Dormand Prince 7, 8  method
af::array AdaptiveRungeKutta::DP78(const State& state, const double dt, double& err_) {
    State tempstate = state;

    double a[14][13] = {{0.}};
    double b[14] = {0.};
    double bhat[14] = {0.};
    double c[14] = {0.};

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

    // Stage 1
    af::array k1 = dt * f(tempstate);

    // Stage 2
    tempstate.t = state.t + c[2] * dt;
    tempstate.m = state.m + a[2][1] * k1;
    af::array k2 = dt * f(tempstate);

    // Stage 3
    tempstate.t = state.t + c[3] * dt;
    tempstate.m = state.m + a[3][1] * k1 + a[3][2] * k2;
    af::array k3 = dt * f(tempstate);

    // Stage 4
    tempstate.t = state.t + c[4] * dt;
    tempstate.m = state.m + a[4][1] * k1 + a[4][2] * k2 + a[4][3] * k3;
    af::array k4 = dt * f(tempstate);

    // Stage 5
    tempstate.t = state.t + c[5] * dt;
    tempstate.m = state.m + a[5][1] * k1 + a[5][2] * k2 + a[5][3] * k3 + a[5][4] * k4;
    af::array k5 = dt * f(tempstate);

    // Stage 6
    tempstate.t = state.t + c[6] * dt;
    tempstate.m = state.m + a[6][1] * k1 + a[6][2] * k2 + a[6][3] * k3 + a[6][4] * k4 + a[6][5] * k5;
    af::array k6 = dt * f(tempstate);

    // Stage 7
    tempstate.t = state.t + c[7] * dt;
    tempstate.m = state.m + a[7][1] * k1 + a[7][2] * k2 + a[7][3] * k3 + a[7][4] * k4 + a[7][5] * k5 + a[7][6] * k6;
    af::array k7 = dt * f(tempstate);

    // Stage 8
    tempstate.t = state.t + c[8] * dt;
    tempstate.m = state.m + a[8][1] * k1 + a[8][2] * k2 + a[8][3] * k3 + a[8][4] * k4 + a[8][5] * k5 + a[8][6] * k6 +
                  a[8][7] * k7;
    af::array k8 = dt * f(tempstate);

    // Stage 9
    tempstate.t = state.t + c[9] * dt;
    tempstate.m = state.m + a[9][1] * k1 + a[9][2] * k2 + a[9][3] * k3 + a[9][4] * k4 + a[9][5] * k5 + a[9][6] * k6 +
                  a[9][7] * k7 + a[9][8] * k8;
    af::array k9 = dt * f(tempstate);

    // Stage 10
    tempstate.t = state.t + c[10] * dt;
    tempstate.m = state.m + a[10][1] * k1 + a[10][2] * k2 + a[10][3] * k3 + a[10][4] * k4 + a[10][5] * k5 +
                  a[10][6] * k6 + a[10][7] * k7 + a[10][8] * k8 + a[10][9] * k9;
    af::array k10 = dt * f(tempstate);

    // Stage 11
    tempstate.t = state.t + c[11] * dt;
    tempstate.m = state.m + a[11][1] * k1 + a[11][2] * k2 + a[11][3] * k3 + a[11][4] * k4 + a[11][5] * k5 +
                  a[11][6] * k6 + a[11][7] * k7 + a[11][8] * k8 + a[11][9] * k9 + a[11][10] * k10;
    af::array k11 = dt * f(tempstate);

    // Stage 12
    tempstate.t = state.t + c[12] * dt;
    tempstate.m = state.m + a[12][1] * k1 + a[12][2] * k2 + a[12][3] * k3 + a[12][4] * k4 + a[12][5] * k5 +
                  a[12][6] * k6 + a[12][7] * k7 + a[12][8] * k8 + a[12][9] * k9 + a[12][10] * k10 + a[12][11] * k11;
    af::array k12 = dt * f(tempstate);

    // Stage 13
    tempstate.t = state.t + c[13] * dt;
    tempstate.m = state.m + a[13][1] * k1 + a[13][2] * k2 + a[13][3] * k3 + a[13][4] * k4 + a[13][5] * k5 +
                  a[13][6] * k6 + a[13][7] * k7 + a[13][8] * k8 + a[13][9] * k9 + a[13][10] * k10 + a[13][11] * k11 +
                  a[13][12] * k12;
    af::array k13 = dt * f(tempstate);

    af::array sumbk = bhat[1] * k1 + bhat[2] * k2 + bhat[3] * k3 + bhat[4] * k4 + bhat[5] * k5 + bhat[6] * k6 +
                      bhat[7] * k7 + bhat[8] * k8 + bhat[9] * k9 + bhat[10] * k10 + bhat[11] * k11 + bhat[12] * k12 +
                      bhat[13] * k13;
    af::array rk_error = sumbk - (b[1] * k1 + b[2] * k2 + b[3] * k3 + b[4] * k4 + b[5] * k5 + b[6] * k6 + b[7] * k7 +
                                  b[8] * k8 + b[9] * k9 + b[10] * k10 + b[11] * k11 + b[12] * k12 + b[13] * k13);
    err_ = max_4d_abs(rk_error / controller_.givescale(max(state.m, state.m + sumbk)));
    return sumbk;
}

// Bogacki-Shampine 2/3rd order  with stepsize control
af::array AdaptiveRungeKutta::BS23(const State& state, const double dt, double& err) {
    State tempstate = state;
    af::array k1;

    if (controller_.get_reject() || normalize_ || step_calls_ == 0) {
        k1 = f(tempstate);
    } else {
        k1 = k_FSAL;
    }

    // stage 2
    tempstate.t = state.t + 1. / 2. * dt;
    tempstate.m = state.m + dt * (1. / 2. * k1);
    ;
    af::array k2 = f(tempstate);

    // stage 3
    tempstate.t = state.t + 3. / 4. * dt;
    tempstate.m = state.m + dt * 3. / 4. * k2;
    ;
    af::array k3 = f(tempstate);

    af::array sumbk = dt * (2. / 9. * k1 + 1. / 3. * k2 + 4. / 9. * k3);

    // stage 4
    tempstate.t = state.t + dt;
    tempstate.m = state.m + sumbk;
    ;
    k_FSAL = f(tempstate);

    af::array rk_error = sumbk - dt * (7. / 24. * k1 + 1. / 4. * k2 + 1. / 3. * k3 + 1. / 8. * k_FSAL);
    err = max_4d_abs(rk_error / controller_.givescale(max(state.m, state.m + sumbk)));
    return sumbk;
}

//// FOR DP and BS, check why error is rising at the beginning of analytical
/// example and then decreases again, maybe use different starting values
//////TODO far too high error in integration test
//// Dormand-Prince 4/5 method
//// FOR DP and BS, check why error is rising at the beginning of analytical
/// example and then decreases again, maybe use different starting values
// af::array AdaptiveRungeKutta::DP45(const State& state, const double dt,
// double& err_)
//{
//    State tempstate=state;
//    // Iterating over a-matrix and calculating k[i]s
//    af::array k[8];
//    const int s = 7;
//
//    double a[8][7]={{0}};
//    double e[8]={0};
//    double c[8]={0};
//
//    c[2]=0.2, c[3]=0.3, c[4]=0.8, c[5]=8.0/9.0, c[6]=1, c[7]=1;
//    e[1]=71.0/57600.0, e[3]=-71.0/16695.0, e[4]=71.0/1920.0,
//    e[5]=-17253.0/339200.0, e[6]=22.0/525.0, e[7]=-1.0/40.0, a[2][1]=0.2,
//    a[3][1]=3.0/40.0, a[3][2]=9.0/40.0,
//    a[4][1]=44.0/45.0, a[4][2]=-56.0/15.0, a[4][3]=32.0/9.0,
//    a[5][1]=19372.0/6561.0, a[5][2]=-25360.0/2187.0, a[5][3]=64448.0/6561.0,
//    a[5][4]=-212.0/729.0, a[6][1]=9017.0/3168.0, a[6][2]=-355.0/33.0,
//    a[6][3]=46732.0/5247.0, a[6][4]=49.0/176.0, a[6][5]=-5103.0/18656.0,
//    a[7][1]=35.0/384.0, a[7][3]=500.0/1113.0, a[7][4]=125.0/192.0,
//    a[7][5]=-2187.0/6784.0, a[7][6]=11.0/84.0;
//
//    // Stage 1
//    if( controller_.get_reject() || normalize_ || step_calls_ == 0)
//    {
//        k[1]   = dt * f(tempstate);
//    }
//    else
//    {
//        k[1]=k_FSAL;
//    }
//    // Stages 2-7
//    for(int i=2;i<=s;i++){
//        af::array rktemp=af::constant(0.0, state.m.dims(0), state.m.dims(1),
//        state.m.dims(2), state.m.dims(3), f64); for(int j=1;j<i;j++){
//            rktemp+=a[i][j] * k[j];
//        }
//        tempstate.t = state.t + c[i];
//        tempstate.m = state.m + rktemp;
//        k[i]= dt * f(tempstate);
//    }
//    //Local extrapolation using 5th order approx
//    af::array sumbk=af::constant(0.0, state.m.dims(0), state.m.dims(1),
//    state.m.dims(2), state.m.dims(3), f64); for(int i=1;i<s;i++){
//        sumbk+=a[s][i]*k[i];
//    }
//    //Error estimation using 4th order approx
//    af::array rk_error=af::constant(0.0, state.m.dims(0), state.m.dims(1),
//    state.m.dims(2), state.m.dims(3), f64); for(int i=1;i<=s;i++){
//        rk_error+=e[i]*k[i];
//    }
//    k_FSAL=k[7];
//    //!!!Note: here e is already the difference between the ususal b and
//    bhat!!!! (no rk_error=sumbk-rk_error)
//    err_=max_4d_abs(rk_error/controller_.givescale(max(state.m, state.m+sumbk)));
//    return sumbk;
//}

//// Bogacki 4, 5 method with sigle error andstepsize control
////TODO far too high error in integration test
// af::array AdaptiveRungeKutta::BS45(const State& state, const double dt ,
// double& err_)
//{
//    State tempstate=state;
//
//    af::array k[9];
//    const int s=8;
//    double a[9][8]={{0.}};
//    double b[9]={0.};
//    double c[9]={0.};
//
//    a[2][1] = 1.0e0/6.0e0;
//    a[3][1] = 2.e0/27.e0;
//    a[3][2] = 4.e0/27.e0;
//    a[4][1] = 183.e0/1372.e0;
//    a[4][2] = -162.e0/343.e0;
//    a[4][3] = 1053.e0/1372.e0;
//    a[5][1] = 68.e0/297.e0;
//    a[5][2] = -4.e0/11.e0;
//    a[5][3] = 42.e0/143.e0;
//    a[5][4] = 1960.e0/3861.e0;
//    a[6][1] = 597.e0/22528.e0;
//    a[6][2] = 81.e0/352.e0;
//    a[6][3] = 63099.e0/585728.e0;
//    a[6][4] = 58653.e0/366080.e0;
//    a[6][5] = 4617.e0/20480.e0;
//    a[7][1] = 174197.e0/959244.e0;
//    a[7][2] = -30942.e0/79937.e0;
//    a[7][3] = 8152137.e0/19744439.e0;
//    a[7][4] = 666106.e0/1039181.e0;
//    a[7][5] = -29421.e0/29068.e0;
//    a[7][6] = 482048.e0/414219.e0;
//    a[8][1] = 587.e0/8064.e0;
//    a[8][2] = 0.e0;
//    a[8][3] = 4440339.e0/15491840.e0;
//    a[8][4] = 24353.e0/124800.e0;
//    a[8][5] = 387.e0/44800.e0;
//    a[8][6] = 2152.e0/5985.e0;
//    a[8][7] = 7267.e0/94080.e0;
//  //C
//  //C
//    c[1] = 0.e0;
//    c[2] = 1.e0/6.e0;
//    c[3] = 2.e0/9.e0;
//    c[4] = 3.e0/7.e0;
//    c[5] = 2.e0/3.e0;
//    c[6] = 3.e0/4.e0;
//    c[7] = 1.e0;
//    c[8] = 1.e0;
//  //C  The coefficients B(*) refer to the formula of order 4.
//    b[1] = 2479.e0/34992.e0;
//    b[2] = 0.e0;
//    b[3] = 123.e0/416.e0;
//    b[4] = 612941.e0/3411720.e0;
//    b[5] = 43.e0/1440.e0;
//    b[6] = 2272.e0/6561.e0;
//    b[7] = 79937.e0/1113912.e0;
//    b[8] = 3293.e0/556956.e0;
//    const bool llg_wasnormalized=true;//TODO
//    if(controller_.get_reject() || (( controller_.get_counter_reject() +
//    controller_.get_counter_accepted()) <=0) || llg_wasnormalized){
//    //if(reject || calls==0 || llg_wasnormalized){
//    //Note: in generalized rkcall use: if(FSAL==false || reject || calls==0 ||
//    llg_wasnormalized){
//      k[1]   =  f(state);
//      //k[1]   =  f(t, m);
//      }
//    else
//      k[1]=k[s];
//    for(int i=2;i<=s;i++){
//        af::array rktemp=af::constant(0.0, state.m.dims(0), state.m.dims(1),
//        state.m.dims(2), state.m.dims(3), f64);
//      for(int j=1;j<i;j++){
//        rktemp+=a[i][j] * k[j];
//      }
//      rktemp*=dt;
//      tempstate.t = state.t + c[i];
//      tempstate.m = state.m + rktemp;
//      k[i]= f(tempstate);
//      //k[i]= f(t + c[i], m + rktemp);
//    }
//
//    af::array sumbk=af::constant(0.0, state.m.dims(0), state.m.dims(1),
//    state.m.dims(2), state.m.dims(3), f64); for(int i=1;i<s;i++){
//      sumbk+=a[s][i]*k[i];
//    }
//    sumbk*=dt;
//
//
//    af::array rk_error=af::constant(0.0, state.m.dims(0), state.m.dims(1),
//    state.m.dims(2), state.m.dims(3), f64); for(int i=1;i<=s;i++){
//      rk_error+=b[i]*k[i];
//    }
//    rk_error*=dt;
//    rk_error=sumbk-rk_error;
//    err_=max_4d_abs(rk_error/controller_.givescale(max(state.m, state.m+sumbk)));
//
//    return sumbk;
//}
} // namespace magnumafcpp
