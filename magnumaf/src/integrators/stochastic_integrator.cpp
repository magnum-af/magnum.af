#include "stochastic_integrator.hpp"
#include "util/util.hpp"
#include <chrono>
#include <utility>


namespace magnumaf {

af::array Stochastic_Integrator::Heun(const State& state) {
    const double D = (alpha * constants::kb * T) /
                     (constants::gamma * constants::mu0 * state.Ms * state.mesh.dx * state.mesh.dy * state.mesh.dz);
    const af::array h_th = sqrt((2. * D) / dt) * af::randn(mesh::dims_v(state.mesh), f64,
                                                           rand_engine); // Random thermal field at t+dt/2
    af::array k1 = dt * stochfdmdt(state, h_th_prev);
    af::array k2 = dt * stochfdmdt(state + k1, h_th);
    h_th_prev = h_th;
    return (k1 + k2) / 2.;
}

af::array Stochastic_Integrator::SemiImplicitHeun(const State& state) {
    const double D = (alpha * constants::kb * T) /
                     (constants::gamma * constants::mu0 * state.Ms * state.mesh.dx * state.mesh.dy * state.mesh.dz);
    const af::array h_th_init =
        sqrt((2. * D) / dt) * randn(mesh::dims_v(state.mesh), f64, rand_engine); // Random thermal field at t
    const af::array h_th = sqrt((2. * D) / dt) * randn(mesh::dims_v(state.mesh), f64,
                                                       rand_engine); // Random thermal field at t+dt/2
    af::array m1 = dt / 2. * stochfdmdt(state, h_th_init);
    af::array m2 = dt / 2. * stochfdmdt(state + m1, h_th);
    af::array m3 = dt / 2. * stochfdmdt(state + m2, h_th);
    af::array m4 = dt / 2. * stochfdmdt(state + m3, h_th);
    af::array m5 = dt / 2. * stochfdmdt(state + m4, h_th);
    return dt * stochfdmdt(state + m5, h_th);
}

af::array Stochastic_Integrator::detRK4(const State& state) const {
    af::array k1 = dt * detfdmdt(state);
    af::array k2 = dt * detfdmdt(state + 1. / 2. * k1);
    af::array k3 = dt * detfdmdt(state + 1. / 2. * k2);
    af::array k4 = dt * detfdmdt(state + k3);
    return (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
}

void Stochastic_Integrator::step(State& state) { // TODO remove dt as parameter here, inconsistency between
                                                 // Heun/SemiHeun
    af::timer timer_stoch = af::timer::start();
    if (mode == 0) {
        state.m += Heun(state);
    } else if (mode == 1) {
        state.m += SemiImplicitHeun(state);
    } else if (mode == 2) {
        state.m += detRK4(state);
    }
    state.m = util::normalize(state.m);
    state.t += dt;
    calls++;
    timer += af::timer::stop(timer_stoch);
    // std::cout<<" TIME  = "<<timer<<std::endl;
}

int get_mode(const std::string& smode) {
    // Setting int mode for usage in void step(...)
    if (smode == "Heun" || smode == "0") {
        return 0;
    } else if (smode == "SemiHeun" || smode == "1") {
        return 1;
    } else if (smode == "detRK4" || smode == "2") {
        return 2;
    } else {
        std::cout << "ERROR: Constructor Stochastic_Ingetrator: Integrator "
                     "Mode not recognized, using Heun (0)"
                  << std::endl;
        return 0;
    }
}

Stochastic_Integrator::Stochastic_Integrator(double alpha, double T, double dt, const State& state,
                                             std::vector<std::unique_ptr<FieldTerm>> fieldterms_in, const std::string& smode)
    : alpha(alpha), T(T), dt(dt), fieldterms(std::move(fieldterms_in)), m_prev(state.m), mode(get_mode(smode)) {
    const double D = (alpha * constants::kb * T) /
                     (constants::gamma * constants::mu0 * state.Ms * state.mesh.dx * state.mesh.dy * state.mesh.dz);
    unsigned long long int seed =
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
            .count();
    rand_engine = af::randomEngine(AF_RANDOM_ENGINE_DEFAULT, seed);
    h_th_prev = sqrt((2. * D) / dt) * randn(mesh::dims_v(state.mesh), f64,
                                            rand_engine); // Initial random thermal field at t=0
}

double Stochastic_Integrator::cpu_time() {
    double cpu_time = 0.;
    for (auto & fieldterm : fieldterms) {
        cpu_time += fieldterm->elapsed_eval_time();
    }
    return cpu_time;
}
} // namespace magnumaf
