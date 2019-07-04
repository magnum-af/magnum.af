#include "atomistic_anisotropy.hpp"
#include "../func.hpp"
//Ref Master Thesis Stifano
//eq (19)
//Han=N*K/2 - K/2 Sum_i(m_i*ez)^2
using namespace af;


double AtomisticUniaxialAnisotropyField::E(const State& state){
  return - (constants::mu0*state.material.p/2.) * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3));
}

double AtomisticUniaxialAnisotropyField::E(const State& state, const af::array& h){
  return - (constants::mu0*state.material.p/2.) * afvalue(sum(sum(sum(sum(h * state.m, 0), 1), 2), 3));
}


AtomisticUniaxialAnisotropyField::AtomisticUniaxialAnisotropyField (const Mesh& mesh, const Material& material){
  //Normal vector
  double norm=sqrt(pow(material.K_atom_axis[0], 2)+ pow(material.K_atom_axis[1], 2) + pow(material.K_atom_axis[2], 2));
  eu=array(mesh.n0, mesh.n1, mesh.n2, 3, f64);
  eu(span, span, span, 0)=material.K_atom_axis[0]/norm;
  eu(span, span, span, 1)=material.K_atom_axis[1]/norm;
  eu(span, span, span, 2)=material.K_atom_axis[2]/norm;
}

array AtomisticUniaxialAnisotropyField::h(const State& state){
  timer_anisotropy = timer::start();
  array anisotropy = eu*state.m;
  anisotropy=sum(anisotropy, 3);
  anisotropy=tile(anisotropy, 1, 1, 1, 3);

  if(state.material.afsync) af::sync();
  cpu_time += timer::stop(timer_anisotropy);
  return  2*state.material.K_atom/(constants::mu0*state.material.p) * anisotropy * eu;
}















////Energy calculation
////Edemag=-mu0/2 integral(M . Hdemag) dx
//double AtomisticUniaxialAnisotropyField::E(const State& state){
//  return - constants::mu0  * state.material.p/2. *  afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3));
//   //return -constants::mu0/2. * state.Ms * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
//}
//
//
////AtomisticUniaxialAnisotropyField::AtomisticUniaxialAnisotropyField (const Mesh& mesh, const Material& material){
////  //Normal vector
////  double norm=sqrt(pow(material.Ku1_axis[0], 2)+ pow(material.Ku1_axis[1], 2) + pow(material.Ku1_axis[2], 2));
////  eu=array(mesh.n0, mesh.n1, mesh.n2, 3, f64);
////  eu(span, span, span, 0)=material.Ku1_axis[0]/norm;
////  eu(span, span, span, 1)=material.Ku1_axis[1]/norm;
////  eu(span, span, span, 2)=material.Ku1_axis[2]/norm;
////}
//
//array AtomisticUniaxialAnisotropyField::h(const State& state){
//  timer_anisotropy = timer::start();
//  double norm=sqrt(pow(state.material.Ku1_axis[0], 2)+ pow(state.material.Ku1_axis[1], 2) + pow(state.material.Ku1_axis[2], 2));
//  array anisotropy=array(state.mesh.n0, state.mesh.n1, state.mesh.n2, 3, f64);
//  //array anisotropy = eu*state.m;
//  anisotropy(span, span, span, 0)=state.material.Ku1_axis[0]/norm * state.m(span, span, span, 0);
//  anisotropy(span, span, span, 1)=state.material.Ku1_axis[1]/norm * state.m(span, span, span, 1);
//  anisotropy(span, span, span, 2)=state.material.Ku1_axis[2]/norm * state.m(span, span, span, 2);
//
//
//  anisotropy*=anisotropy;
//  anisotropy=sum(anisotropy, 3);
//  anisotropy=tile(anisotropy, 1, 1, 1, 3);
//
//  if(state.material.afsync) sync();
//  cpu_time += timer::stop(timer_anisotropy);
//  return  state.mesh.n0*state.mesh.n1*state.mesh.n2*3*state.material.Ku1/2. - state.material.Ku1/2. * anisotropy;
//}
