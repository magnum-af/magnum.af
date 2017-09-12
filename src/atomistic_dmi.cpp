#include "atomistic_dmi.hpp"
using namespace af;
//Energy calculation
//E=-mu0/2 integral(M . H) dx
double ATOMISTIC_DMI::E(const State& state){
  return - state.param.mu0/2. * state.param.p * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)); 
}

ATOMISTIC_DMI::ATOMISTIC_DMI (const Mesh& mesh, const Param& param){
  //initialize finite difference first order derivative filter
  double norm=sqrt(pow(param.D_atom_axis[0],2)+ pow(param.D_atom_axis[1],2) + pow(param.D_atom_axis[2], 2));
  n=array(mesh.n0,mesh.n1,mesh.n2,3,f64);
  n(span,span,span,0)=param.D_atom_axis[0]/norm;
  n(span,span,span,1)=param.D_atom_axis[1]/norm;
  n(span,span,span,2)=param.D_atom_axis[2]/norm;
  //print("n",n);

  filtr_fd1=constant(0.0,3,3,3,3,f64);
  //dmx/dx
  filtr_fd1(0,1,1,0)= 1;
  filtr_fd1(2,1,1,0)=-1;

  //dmy/dy
  filtr_fd1(1,0,1,1)= 1;
  filtr_fd1(1,2,1,1)=-1;

  //dmz/dz
  filtr_fd1(1,1,0,2)= 1;
  filtr_fd1(1,1,2,2)=-1;
}

array ATOMISTIC_DMI::h(const State& state){
  timer_dmi = timer::start();
  array first = convolve(state.m,filtr_fd1,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
  //Make Divergence
  first=sum(first,3);
  first=tile(first,1,1,1,3);
  first=n*first;

  //Dot product
  array second = sum(n*state.m,3);
  //Expand for fd1 convolution
  second = tile(second,1,1,1,3);
  second = convolve(second,filtr_fd1,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
  if(state.param.afsync) sync();
  cpu_time += timer::stop(timer_dmi);
  return -state.param.D_atom/(2.*state.param.mu0*state.param.p) * (first-second);
  //return -state.param.D/2. * res;
}





////VERSION FOR fixed e_z=(0,0,1)
//#include "atomistic_dmi.hpp"
//using namespace af;
////Energy calculation
////E=-mu0/2 integral(M . H) dx
//double ATOMISTIC_DMI::E(const State& state){
//  return - state.param.mu0/2. * state.param.p * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)); 
//  //return -state.param.mu0/2. * state.param.p * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)); 
//  //return -state.param.mu0/2. * state.param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz; 
//}
//
//ATOMISTIC_DMI::ATOMISTIC_DMI (const Mesh& mesh, const Param& param){
//  //initialize finite difference first order derivative filter
//  double norm=sqrt(pow(param.D_axis[0],2)+ pow(param.D_axis[1],2) + pow(param.D_axis[2], 2));
//  n=array(mesh.n0,mesh.n1,mesh.n2,3,f64);
//  n(span,span,span,0)=param.D_axis[0]/norm;
//  n(span,span,span,1)=param.D_axis[1]/norm;
//  n(span,span,span,2)=param.D_axis[2]/norm;
//  //print("n",n);
//
//  //initialize finite difference first order derivative filter
//  filtr_atom_dmi_x=constant(0.0,3,3,3,3,f64);
//  filtr_atom_dmi_y=constant(0.0,3,3,3,3,f64);
//  filtr_atom_dmi_z=constant(0.0,3,3,3,3,f64);
//  //dmx/dx
//  filtr_atom_dmi_x(0,1,1,2)= 1;
//  filtr_atom_dmi_x(2,1,1,2)=-1;
//  //print("filtr",filtr_atom_dmi_x);
//
//  //dmy/dy
//  filtr_atom_dmi_y(1,0,1,2)= 1;
//  filtr_atom_dmi_y(1,2,1,2)=-1;
//  //print("filtr_y",filtr_atom_dmi_y);
//
//  //dmz/dz
//  filtr_atom_dmi_z(0,1,1,0)=-1;
//  filtr_atom_dmi_z(2,1,1,0)= 1;
//  filtr_atom_dmi_z(1,0,1,1)=-1;
//  filtr_atom_dmi_z(1,2,1,1)= 1;
//  //print("filtr_z",filtr_atom_dmi_z);
//}
//
//array ATOMISTIC_DMI::h(const State& state){
//  timer_dmi = timer::start();
//
//  //print("m",state.m);
//  array convx = convolve(state.m,filtr_atom_dmi_x,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
//  //print("convx",convx);
//  array convy = convolve(state.m,filtr_atom_dmi_y,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
//  //print("convy",convy);
//  array convz = convolve(state.m,filtr_atom_dmi_z,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
//  //print("convz",convz);
//  array res = join(3, convx(span,span,span,2), convy(span,span,span,2),convz(span,span,span,0)+convz(span,span,span,1));
//  //print("res",res);
//  if(state.param.afsync) sync();
//  cpu_time += timer::stop(timer_dmi);
//  return -state.param.D_atom/(2.*state.param.mu0*state.param.p) * res;
//  //return -state.param.D/2. * res;
//}
