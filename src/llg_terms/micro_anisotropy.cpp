//Ref Diss Abert p 15 sec 2.3 eq 2.3.28
//Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m) 
//With higher order (not implemented): Hu(r)=2 Ku1 /(mu0 Ms) eu ( eu . m) ( + 4 Ku2 /(mu0 Ms) eu ( eu . m)^3
#include "micro_anisotropy.hpp"

//Energy calculation
//Edemag=-mu0/2 integral(M . Hdemag) dx
double ANISOTROPY::E(const State& state){
  return -param.mu0/2. * param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz; 
}

double ANISOTROPY::E(const State& state, const af::array& h){
  return -param.mu0/2. * param.ms * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz; 
}


ANISOTROPY::ANISOTROPY (Mesh meshin, Param paramin) : param(paramin),mesh(meshin){
  //Normal vector
  double norm=sqrt(pow(param.Ku1_axis[0],2)+ pow(param.Ku1_axis[1],2) + pow(param.Ku1_axis[2], 2));
  eu = af::array(mesh.n0,mesh.n1,mesh.n2,3,f64);
  eu(af::span,af::span,af::span,0)=param.Ku1_axis[0]/norm;
  eu(af::span,af::span,af::span,1)=param.Ku1_axis[1]/norm;
  eu(af::span,af::span,af::span,2)=param.Ku1_axis[2]/norm;
}

af::array ANISOTROPY::h(const State& state){
  timer_anisotropy = af::timer::start();
  af::array anisotropy = eu*state.m;
  anisotropy=af::sum(anisotropy,3);
  anisotropy=af::tile(anisotropy,1,1,1,3);

  if(param.afsync) af::sync();
  cpu_time += af::timer::stop(timer_anisotropy);
  return  2.* param.Ku1/(param.mu0 * param.ms) * (eu* anisotropy);
}










//  print("eu",eu);
//  array m = constant(0.0,mesh.n0,mesh.n1,mesh.n2,3,f64);
//  m(seq(1,end-1),span,span,0) = constant(1.0,mesh.n0-2,mesh.n1,mesh.n2,1,f64);
//  m(0,span,span,1 ) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
//  m(-1,span,span,1) = constant(1.0,1,mesh.n1,mesh.n2,1,f64);
//  array anisotropy = eu*m;
//  print("a",anisotropy);
//  anisotropy=sum(anisotropy,3);
//  print("a",anisotropy);
//  anisotropy=tile(anisotropy,1,1,1,3);
//  print("a",anisotropy);
//
