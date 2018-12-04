//Ref Master Thesis Stifano
//eq (19) 
//Han=N*K/2 - K/2 Sum_i(m_i*ez)^2
#include "atomistic_anisotropy.hpp"
using namespace af;


double ATOMISTIC_ANISOTROPY::E(const State& state){
  return - (state.param.mu0*state.param.p/2.) * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)); 
}

double ATOMISTIC_ANISOTROPY::E(const State& state, const af::array& h){
  return - (state.param.mu0*state.param.p/2.) * afvalue(sum(sum(sum(sum(h * state.m,0),1),2),3));
}


ATOMISTIC_ANISOTROPY::ATOMISTIC_ANISOTROPY (const Mesh& mesh, const Param& param){
  //Normal vector
  double norm=sqrt(pow(param.K_atom_axis[0],2)+ pow(param.K_atom_axis[1],2) + pow(param.K_atom_axis[2], 2));
  eu=array(mesh.n0,mesh.n1,mesh.n2,3,f64);
  eu(span,span,span,0)=param.K_atom_axis[0]/norm;
  eu(span,span,span,1)=param.K_atom_axis[1]/norm;
  eu(span,span,span,2)=param.K_atom_axis[2]/norm;
}

array ATOMISTIC_ANISOTROPY::h(const State& state){
  timer_anisotropy = timer::start();
  array anisotropy = eu*state.m;
  anisotropy=sum(anisotropy,3);
  anisotropy=tile(anisotropy,1,1,1,3);

  if(state.param.afsync) sync();
  cpu_time += timer::stop(timer_anisotropy);
  return  2*state.param.K_atom/(state.param.mu0*state.param.p) * anisotropy * eu;
}















////Energy calculation
////Edemag=-mu0/2 integral(M . Hdemag) dx
//double ATOMISTIC_ANISOTROPY::E(const State& state){
//  return - state.param.mu0  * state.param.p/2. *  afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)); 
//   //return -state.param.mu0/2. * state.param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz; 
//}
//
//
////ATOMISTIC_ANISOTROPY::ATOMISTIC_ANISOTROPY (const Mesh& mesh, const Param& param){
////  //Normal vector
////  double norm=sqrt(pow(param.Ku1_axis[0],2)+ pow(param.Ku1_axis[1],2) + pow(param.Ku1_axis[2], 2));
////  eu=array(mesh.n0,mesh.n1,mesh.n2,3,f64);
////  eu(span,span,span,0)=param.Ku1_axis[0]/norm;
////  eu(span,span,span,1)=param.Ku1_axis[1]/norm;
////  eu(span,span,span,2)=param.Ku1_axis[2]/norm;
////}
//
//array ATOMISTIC_ANISOTROPY::h(const State& state){
//  timer_anisotropy = timer::start();
//  double norm=sqrt(pow(state.param.Ku1_axis[0],2)+ pow(state.param.Ku1_axis[1],2) + pow(state.param.Ku1_axis[2], 2));
//  array anisotropy=array(state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
//  //array anisotropy = eu*state.m;
//  anisotropy(span,span,span,0)=state.param.Ku1_axis[0]/norm * state.m(span,span,span,0);
//  anisotropy(span,span,span,1)=state.param.Ku1_axis[1]/norm * state.m(span,span,span,1);
//  anisotropy(span,span,span,2)=state.param.Ku1_axis[2]/norm * state.m(span,span,span,2);
//  
//
//  anisotropy*=anisotropy;
//  anisotropy=sum(anisotropy,3);
//  anisotropy=tile(anisotropy,1,1,1,3);
//
//  if(state.param.afsync) sync();
//  cpu_time += timer::stop(timer_anisotropy);
//  return  state.mesh.n0*state.mesh.n1*state.mesh.n2*3*state.param.Ku1/2. - state.param.Ku1/2. * anisotropy;
//}
