#include "stochastic_integrator.hpp"
using namespace af;

af::array Stochastic_Integrator::Heun(const State& state, const double dt)
{
    const double D = (material.alpha * constants::kb * material.T)/ (constants::gamma * constants::mu0 * this->Ms * mesh.V);
    const array h_th = sqrt ((2. * D)/dt) * randn(mesh.dims, f64, rand_engine);// Random thermal field at t+dt/2
    af::array k1 = dt * stochfdmdt(state, h_th_prev);
    af::array k2 = dt * stochfdmdt(state + k1, h_th);
    h_th_prev = h_th;
    return (k1 + k2 )/2.;
}

af::array Stochastic_Integrator::SemiImplicitHeun(const State& state, const double dt)
{
    const double D = (material.alpha * constants::kb * material.T)/ (constants::gamma * constants::mu0 * this->Ms * mesh.V);
    const array h_th_init = sqrt ((2. * D)/dt) * randn(mesh.dims, f64, rand_engine);// Random thermal field at t
    const array h_th = sqrt ((2. * D)/dt) * randn(mesh.dims, f64, rand_engine);// Random thermal field at t+dt/2
    af::array m1= dt/2. * stochfdmdt (state   , h_th_init);
    af::array m2= dt/2. * stochfdmdt (state + m1, h_th);
    af::array m3= dt/2. * stochfdmdt (state + m2, h_th);
    af::array m4= dt/2. * stochfdmdt (state + m3, h_th);
    af::array m5= dt/2. * stochfdmdt (state + m4, h_th);
    return   dt * stochfdmdt (state + m5, h_th);
}

af::array Stochastic_Integrator::detRK4(const State& state, const double dt)
{
    af::array k1   =  dt * detfdmdt(state                               );
    af::array k2   =  dt * detfdmdt(state + 1./2.*k1                    );
    af::array k3   =  dt * detfdmdt(state            + 1./2.*k2         );
    af::array k4   =  dt * detfdmdt(state                       +    k3 );
    return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
}

void Stochastic_Integrator::step(State& state, const double dt){//TODO remove dt as parameter here, inconsistency between Heun/SemiHeun
    timer_stoch = timer::start();
    if (mode == 0){ 
        state.m += Heun(state,dt);
    }
    else if (mode == 1){ 
        state.m += SemiImplicitHeun(state,dt);
    }
    else if (mode == 2){
        state.m += detRK4(state,dt);
    }
    state.m = renormalize(state.m);
    state.t+=dt;
    calls ++;
    time += timer::stop(timer_stoch);
    //std::cout<<" TIME  = "<<time<<std::endl;
}

Stochastic_Integrator::Stochastic_Integrator (State state, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in, const double dt, std::string smode):
  Fieldterms(Fieldterms_in),  material(state.material), mesh(state.mesh), Ms(state.Ms), m_prev(state.m)
{
    const double D = (material.alpha * constants::kb * material.T)/ (constants::gamma * constants::mu0 * state.Ms * mesh.V);
    unsigned long long int seed = std::chrono::duration_cast< std::chrono::nanoseconds >( std::chrono::system_clock::now().time_since_epoch()).count();
    rand_engine=af::randomEngine(af::randomEngine(AF_RANDOM_ENGINE_DEFAULT, seed));
    h_th_prev = sqrt ((2. * D)/dt) * randn(mesh.dims, f64, rand_engine);// Initial random thermal field at t=0

    //Setting int mode for usage in void step(...)
    if (smode == "Heun" || smode == "0"){ mode = 0;}
    else if (smode == "SemiHeun" || smode == "1"){ mode = 1;}
    else if (smode == "detRK4" || smode == "2"){ mode = 2;}
    else {std::cout<< "ERROR: Constructor Stochastic_Ingetrator: Integrator Mode not recognized"<<std::endl;}
}

double Stochastic_Integrator::cpu_time(){
    double cpu_time = 0.;
    for(unsigned i=0;i<Fieldterms.size();++i){
        cpu_time+=Fieldterms[i]->get_cpu_time();
    }
    return cpu_time;
}
