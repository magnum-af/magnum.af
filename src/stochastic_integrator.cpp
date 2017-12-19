#include "stochastic_integrator.hpp"
Stochastic_Integrator::Stochastic_Integrator (State in, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in, const double dt) : Fieldterms(Fieldterms_in),  param(in.param), mesh(in.mesh), m_prev(in.m){
    const double D = (param.alpha * param.kb * param.T)/ (param.gamma * param.mu0 * param.ms * mesh.V);
    h_th_prev = sqrt ((2. * D)/dt) * randn(mesh.dims, f64);// Initial random thermal field at t=0
}

double Stochastic_Integrator::cpu_time(){
  double cpu_time = 0.;
  for(unsigned i=0;i<Fieldterms.size();++i){
    cpu_time+=Fieldterms[i]->get_cpu_time();
  }
  return cpu_time;
}

template <class T>  T Stochastic_Integrator::StochSemiImplicitHeun(const T& m, const double dt)
{
    const double D = (param.alpha * param.kb * param.T)/ (param.gamma * param.mu0 * param.ms * mesh.V);
    const array h_th_init = sqrt ((2. * D)/dt) * randn(mesh.dims, f64);// Random thermal field at t
    const array h_th = sqrt ((2. * D)/dt) * randn(mesh.dims, f64);// Random thermal field at t+dt/2
    T m1= dt/2. * stochfdmdt (m   , h_th_init);
    T m2= dt/2. * stochfdmdt (m+m1, h_th);
    T m3= dt/2. * stochfdmdt (m+m2, h_th);
    T m4= dt/2. * stochfdmdt (m+m3, h_th);
    T m5= dt/2. * stochfdmdt (m+m4, h_th);
    return   dt * stochfdmdt (m+m5, h_th);
}

template <class T>  T Stochastic_Integrator::StochHeun(const T& m, const double dt)
{
    const double D = (param.alpha * param.kb * param.T)/ (param.gamma * param.mu0 * param.ms * mesh.V);
    const array h_th = sqrt ((2. * D)/dt) * randn(mesh.dims, f64);// Random thermal field at t+dt/2
    T k1 = dt * stochfdmdt(m, h_th_prev);
    T k2 = dt * stochfdmdt(m + k1, h_th);
    h_th_prev = h_th;
    return (k1 + k2 )/2.;
}

template <class T>  T Stochastic_Integrator::SemiImplicitHeun(const T& m,  const double dt)
{
    T m1= dt/2. * fdmdt (m);
    T m2= dt/2. * fdmdt (m+m1);
    T m3= dt/2. * fdmdt (m+m2);
    T m4= dt/2. * fdmdt (m+m3);
    T m5= dt/2. * fdmdt (m+m4);
    return  fdmdt( m + m5) * dt;
}

//array rk4(array (*fdmdt)(const array&)  ,const array& m, const double dt)
template <class T>  T Stochastic_Integrator::rk4(const T& m, const double dt)
{
    T k1   =  dt * fdmdt(m                               );
    T k2   =  dt * fdmdt(m + 1./2.*k1                    );
    T k3   =  dt * fdmdt(m            + 1./2.*k2         );
    T k4   =  dt * fdmdt(m                       +    k3 );
    return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
}

void Stochastic_Integrator::step(State& state, const double dt){
    timer_stoch = timer::start();
    //state.m += rk4(state.m,dt);
    //state.m += SemiImplicitHeun(state.m,dt);
    //state.m += StochSemiImplicitHeun(state.m,dt);
    state.m += StochHeun(state.m,dt);
    state.m = renormalize(state.m);
    state.t+=dt;
    calls ++;
    time += timer::stop(timer_stoch);
    //std::cout<<" TIME  = "<<time<<std::endl;
}


//template <class T>  T Stochastic_Integrator::StochSemiImplicitHeun(const T& m, const double dt)
//{
//    T m_it = 1.5 * m - m_prev/2.;
//    m_prev = m;
//    const double D = (param.alpha * param.kb * param.T)/ (param.gamma * param.mu0 * param.ms * mesh.V);
//    const array h_th = sqrt ((2. * D)/dt) * randn(mesh.dims, f64);
//    for (int i = 0; i < 5; i++){
//        m_it = m + stochfdmdt(m_it, h_th) * dt/2.;
//    }
//    return  stochfdmdt(m_it, h_th) * dt;
//}


//template <class T>  T Stochastic_Integrator::SemiImplicitHeun(const T& m, const double dt)
//{
//    T m_it = 1.5 * m - m_prev/2.;
//    m_prev = m;
//    for (int i = 0; i < 5; i++){
//        m_it = m + fdmdt(m_it) * dt/2.;
//    }
//    return  fdmdt(m_it) * dt;
//}

//template <class T>  T Stochastic_Integrator::SemiImplicitHeun(const T& m, const double dt)
//{
//    //T m0= 1.5* m - 0.5 * m_prev;
//    //m_prev = m;
//    T k1= dt/2. * fdmdt (m);
//    T k2= dt/2. * fdmdt (m+k1);
//    T k3= dt/2. * fdmdt (m+k2);
//    return  dt * fdmdt(m+k3);
//    //T k1   =  dt * fdmdt(m                               );
//    //T k2   =  dt * fdmdt(m + 1./2.*k1                    );
//    //T k3   =  dt * fdmdt(m            + 1./2.*k2         );
//    //T k4   =  dt * fdmdt(m                       +    k3 );
//    //return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
//}
//template <class T>  T Stochastic_Integrator::SemiImplicitHeun(const T& m, const double dt)
//{
//    T m0= 1.5* m - 0.5 * m_prev;
//    m_prev = m;
//    T m1= m + dt/2. * fdmdt (m0);
//    T m2= m + dt/2. * fdmdt (m1);
//    T m3= m + dt/2. * fdmdt (m2);
//    T m4= m + dt/2. * fdmdt (m3);
//    T m5= m + dt/2. * fdmdt (m4);
//    return m + dt * fdmdt(m5);
//}

////array rk4(array (*fdmdt)(const array&)  ,const array& m, const double dt)
//array Stochastic_Integrator::rk4(const array& m, const double dt)
//{
//    array k1   =  dt * fdmdt(m                               );
//    array k2   =  dt * fdmdt(m + 1./2.*k1                    );
//    array k3   =  dt * fdmdt(m            + 1./2.*k2         );
//    array k4   =  dt * fdmdt(m                       +    k3 );
//    return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
//}


//array Stochastic_Integrator::StemiImplicitHeun(const array& m, const double dt){
//    timer_stoch = timer::start();
//    array m_it = m + ( m - m_prev)/2.;
//    for (int i = 0; i < 5; i++){
//       m_it = m + fdmdt(m_it) * dt/2.;
//    }
//    m_prev = m;
//    std::cout<<" CPU TIME  = "<<cpu_time()<<"\n"<<std::endl;
//    return m + fdmdt(m_it) * dt;
//}


//array Stochastic_Integrator::fdmdt(const array& m, const array& heff){
//    array cross_temp = cross4(m, heff);
//    return - param.gamma/(1.+pow(param.alpha,2)) * cross4(m, heff) - param.alpha*param.gamma/(1.+pow(param.alpha,2)) * cross4(m, cross4(m, heff));
//}
//
//array Stochastic_Integrator::rk4(const array& m, const double dt)
//{
//  array heff =       fheff(m                                              );
//  array k1   =  dt * fdmdt(m                                        , heff);
//  heff =       fheff(m + 1./2.*k1                                   );
//  array k2   =  dt * fdmdt(m + 1./2.*k1                             , heff);
//  heff =       fheff(m            + 1./2.*k2                        );
//  array k3   =  dt * fdmdt(m            + 1./2.*k2                  , heff);
//  heff =       fheff(m                       +    k3                );
//  array k4   =  dt * fdmdt(m                       +    k3          , heff);
//
//  return            (          k1 +    2.*k2 + 2.*k3 + k4    ) / 6.;
//}

//void Stochastic_Integrator::SemiImplicitHeun(array& m, const double dt){
//array Stochastic_Integrator::StemiImplicitHeun(const array& m, const double dt){
    //Now in constructor
    //if(calls == 0) { 
    //    m_prev = m; 
    //    std::cout<<"CALLS == 0"<<std::endl;
    //}
    //TODO this is somehow super slow
    //
//    array m_it = m + ( m - m_prev)/2.;
//    for (int i = 0; i < 5; i++){
//       m_it = m + fdmdt(m_it) * dt/2.;
//    }
//    m_prev = m;
//    std::cout<<" CPU TIME  = "<<cpu_time()<<"\n"<<std::endl;
//    //return m + fdmdt(m_it) * dt;
//    return  fdmdt(m_it) * dt;
    
//    array m0= 1.5* m - 0.5 * m_prev;
//    m_prev = m;
//    array m1= m + dt/2. * fdmdt (m0);
//    array m2= m + dt/2. * fdmdt (m1);
//    array m3= m + dt/2. * fdmdt (m2);
//    array m4= m + dt/2. * fdmdt (m3);
//    array m5= m + dt/2. * fdmdt (m4);
//    array result = dt * fdmdt(m5);
//    m += result;

    //array result = dt * fdmdt(m5);
    //return result;
    //return dt *fdmdt(m5);
    //return dt *fdmdt(m5);//TODO extremely slow for first ~100 runs, then fast

//    array k1   =  dt * fdmdt(m                               );
//    array k2   =  dt * fdmdt(m + 1./2.*k1                    );
//    array k3   =  dt * fdmdt(m            + 1./2.*k2         );
//    array k4   =  dt * fdmdt(m                       +    k3 );
//  
//    time += timer::stop(timer_stoch);
//    return (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
//}

//array Stochastic_Integrator::SemiHeun(const array& m, const double dt){
    //array m0= 1.5* m - 0.5 * m_prev;
    //m_prev = m;
    //return randu(100,25,1,3,f64);

    //TODO THIS IS SUPER SLOW 
    //return m;

//}