#include "micro_dmi.hpp"
using namespace af;
void showdims(const array& a){
  std::cout<<"Exchange matrix: dims="<<a.dims(0)<<"\t"<<a.dims(1)<<"\t"<<a.dims(2)<<"\t"<<a.dims(3)<<std::endl;
}
//Energy calculation
//E=-mu0/2 integral(M . H) dx
double DMI::E(const State& state){
  return -param.mu0/2. * param.ms * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3)) * mesh.dx * mesh.dy * mesh.dz; 
}

DMI::DMI (Mesh meshin, Param paramin) : param(paramin),mesh(meshin){
  //Normal vector
  double norm=sqrt(pow(param.D_axis[0],2)+ pow(param.D_axis[1],2) + pow(param.D_axis[2], 2));
  n=array(mesh.n0,mesh.n1,mesh.n2,3,f64);
  n(span,span,span,0)=param.D_axis[0]/norm;
  n(span,span,span,1)=param.D_axis[1]/norm;
  n(span,span,span,2)=param.D_axis[2]/norm;
  //print("n",n);

  //initialize finite difference first order derivative filter
  //Central finite difference in 2nd order
  //https://en.wikipedia.org/wiki/Finite_difference_coefficient
  filtr_fd1=constant(0.0,3,3,3,3,f64);
  //dmx/dx
  filtr_fd1(0,1,1,0)= 1 / (2.*mesh.dx);
  filtr_fd1(2,1,1,0)=-1 / (2.*mesh.dx);

  //dmy/dy
  filtr_fd1(1,0,1,1)= 1 / (2.*mesh.dy);
  filtr_fd1(1,2,1,1)=-1 / (2.*mesh.dy);

  //dmz/dz
  filtr_fd1(1,1,0,2)= 1 / (2.*mesh.dz);
  filtr_fd1(1,1,2,2)=-1 / (2.*mesh.dz);
}


array DMI::h(const State& state){
  timer_dmi = timer::start();

  //First: n(div m)
  //Gradient and edges
  array first = convolve(state.m,filtr_fd1,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
  correct_edges(first,state.m);
  //Make Divergence
  first=sum(first,3);
  first=tile(first,1,1,1,3);
  first=n*first;

  //Second: fd1(n . m)
  //Dot product
  array second = sum(n*state.m,3);
  //Expand for fd1 convolution
  second = tile(second,1,1,1,3);
  second = convolve(second,filtr_fd1,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
  correct_edges(second,state.m);

  if(param.afsync) sync();
  cpu_time += timer::stop(timer_dmi);
  return 2.* param.D/(param.mu0*param.ms) * (first-second);//Note: Js=mu0*Ms
}


//EXPERIMENTAL
//This version assumes the same boundaries as in the exchange filed of 70Lines of Numpy
//On the edges, we assume (not existing) outer cells with the same values as the boundary cells
// in the middle ...|-1/2|  0 |1/2|... 
// on the edges         ||  0 |1/2|...
// we want         |-1/2||  0 |1/2|... 
// thus we take         ||-1/2|1/2|...
// so after the convolution we have to add the edges with -1/(2*dx)

void DMI::correct_edges(array& out, const array& in){
  //Lower x edge:
  out( 0,span,span,0)+= -0.5* in( 0,span,span,0)/mesh.dx;
  //Upper x edge:
  out(-1,span,span,0)-= -0.5* in(-1,span,span,0)/mesh.dx;

  //Lower y edge:
  out(span, 0,span,1)+= -0.5* in(span, 0,span,1)/mesh.dy;
  //Upper y edge:
  out(span,-1,span,1)-= -0.5* in(span,-1,span,1)/mesh.dy;

  //z
  if(in.dims(2)==1){
    out(span,span,span,2)=0.;
  }
  else{
    //Lower z edge:
    out(span,span, 0,2)+= -0.5* in(span,span, 0,2)/mesh.dz;
    //Upper z edge:
    out(span,span,-1,2)-= -0.5* in(span,span,-1,2)/mesh.dz;
  }
}


//void DMI::correct_edges(array& out, const array& in){
//  //Lower x bound: after convolve it is:  1/2 * (1)  
//  //i.e.                  filtr_grad(0,1,1,0) * in(0,span,span,0)
//  //We want it to be:                  -1*(0) + 1* (1)
//  //So we take                        1/2*(1) - 1* (0)
//  
//  //Lower x edge:
//  out( 0,span,span,0)+= in( 1,span,span,0)/(2.*mesh.dx) - in( 0,span,span,0)/mesh.dx;
//  //Upper x edge:
//  out(-1,span,span,0)-= in(-2,span,span,0)/(2.*mesh.dx) - in(-1,span,span,0)/mesh.dx;
//
//  //Lower y edge:
//  out(span, 0,span,1)+= in(span, 1,span,1)/(2.*mesh.dy) - in(span, 0,span,1)/mesh.dy;
//  //Upper y edge:                                                                  
//  out(span,-1,span,1)-= in(span,-2,span,1)/(2.*mesh.dy) - in(span,-1,span,1)/mesh.dy;
//
//  //z
//  if(in.dims(2)==1){
//    out(span,span,span,2)=0.;
//  }
//  else{
//    //Lower z edge:
//    out(span,span, 0,2)+= in(span,span, 1,2)/(2.*mesh.dz) - in(span,span, 0,2)/mesh.dz;
//    //Upper z edge:                                                                  
//    out(span,span,-1,2)-= in(span,span,-2,2)/(2.*mesh.dz) - in(span,span,-1,2)/mesh.dz;
//  }
//}















//DMI::DMI (Mesh meshin, Param paramin) : param(paramin),mesh(meshin){
//  //Normal vector
//  n=array(mesh.n0,mesh.n1,mesh.n2,3,f64);
//  n(span,span,span,0)=nx;
//  n(span,span,span,1)=ny;
//  n(span,span,span,2)=nz;
//  //print("n",n);
//
//  //TODO:
//  //initialize filters
//  filtr=constant(0.0,3,3,3,f64);
//  //filtr(1,1,1)= -6 / (pow(mesh.dx,2)+pow(mesh.dy,2)+pow(mesh.dz,2));
//
//  filtr(0,1,1)=-1 / (2.*mesh.dx);
//  filtr(2,1,1)= 1 / (2.*mesh.dx);
//
//  filtr(1,0,1)=-1 / (2.*mesh.dy);
//  filtr(1,2,1)= 1 / (2.*mesh.dy);
//
//  filtr(1,1,0)=-1 / (2.*mesh.dz);
//  filtr(1,1,2)= 1 / (2.*mesh.dz);
//}
//array DMI::h(array m){
//  timer_dmi = timer::start();
//  array first = convolve(m,filtr,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
//  showdims(first);
//  first=n*tile(sum(first,3),1,1,1,3);
//  showdims(first);
//  array second = convolve(tile(sum(n*m,3),1,1,1,3),filtr,AF_CONV_DEFAULT,AF_CONV_SPATIAL);
//  showdims(second);
//
//  if(param.afsync) sync();
//  cpu_time += timer::stop(timer_dmi);
//  return  -2.* param.D/param.Js * (first-second);//TODO Js not set
//}
