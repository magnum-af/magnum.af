#include "atomistic_dmi.hpp"
#include "../func.hpp"

namespace magnumaf{

using namespace af;
//Energy calculation
//E=-mu0/2 integral(M . H) dx
float AtomisticDmiField::E(const State& state){
  return - constants::mu0/2. * state.material.p * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3));
}

float AtomisticDmiField::E(const State& state, const af::array& h){
  return - constants::mu0/2. * state.material.p * afvalue(sum(sum(sum(sum(h * state.m, 0), 1), 2), 3));
}

AtomisticDmiField::AtomisticDmiField (const Mesh& mesh, const Material& material){
    //initialize finite difference first order derivative filter
    float norm=sqrt(pow(material.D_atom_axis[0], 2)+ pow(material.D_atom_axis[1], 2) + pow(material.D_atom_axis[2], 2));
    n=array(mesh.n0, mesh.n1, mesh.n2, 3, f32);
    n(span, span, span, 0)=material.D_atom_axis[0]/norm;
    n(span, span, span, 1)=material.D_atom_axis[1]/norm;
    n(span, span, span, 2)=material.D_atom_axis[2]/norm;
    //print("n", n);

    filtr_fd1=constant(0.0, 3, 3, 3, 3, f32);
    //dmx/dx
    filtr_fd1(0, 1, 1, 0)= 1;
    filtr_fd1(2, 1, 1, 0)=-1;

    //dmy/dy
    filtr_fd1(1, 0, 1, 1)= 1;
    filtr_fd1(1, 2, 1, 1)=-1;

    //dmz/dz
    filtr_fd1(1, 1, 0, 2)= 1;
    filtr_fd1(1, 1, 2, 2)=-1;
}

array AtomisticDmiField::h(const State& state){
    timer_dmi = timer::start();
    array first = convolve(state.m, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    //Make Divergence
    first=sum(first, 3);
    first=tile(first, 1, 1, 1, 3);
    first=n*first;

    //Dot product
    array second = sum(n*state.m, 3);
    //Expand for fd1 convolution
    second = tile(second, 1, 1, 1, 3);
    second = convolve(second, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    if(state.material.afsync) af::sync();
    cpu_time += timer::stop(timer_dmi);
    return state.material.D_atom/(constants::mu0*state.material.p) * (first-second);
    //return -state.material.D/2. * res;
}


////Version with DmiField d pointing only in nz direction
//array AtomisticDmiField::h(const State& state){
//    timer_dmi = timer::start();
//
//    array filterHdx=constant(0., 3, 1, 1, 3, f32);
//    filterHdx(2, 0, 0, 2)=-1;
//    filterHdx(0, 0, 0, 2)= 1;
//    array Hdx = convolve( state.m, filterHdx, AF_CONV_DEFAULT, AF_CONV_SPATIAL)(span, span, span, 2);
//    //print ("filterHdx", filterHdx);
//    //print ("Hdx", Hdx);
//
//    array filterHdy=constant(0., 1, 3, 1, 3, f32);
//    filterHdy(0, 2, 0, 2)=-1;
//    filterHdy(0, 0, 0, 2)= 1;
//    array Hdy = convolve( state.m, filterHdy, AF_CONV_DEFAULT, AF_CONV_SPATIAL)(span, span, span, 2);
//    //print ("Hdy", Hdy);
//
//    array filterHdz=constant(0., 3, 3, 1, 3, f32);
//    filterHdz(2, 1, 0, 0)= 1;
//    filterHdz(0, 1, 0, 0)=-1;
//    filterHdz(1, 2, 0, 1)= 1;
//    filterHdz(1, 0, 0, 1)=-1;
//    array Hdz = convolve( state.m, filterHdz, AF_CONV_DEFAULT, AF_CONV_SPATIAL)(span, span, span, seq(0, 1));
//    Hdz=Hdz(span, span, span, 0)+Hdz(span, span, span, 1);
//    //print ("Hdz", Hdz);
//
//    if(state.material.afsync) sync();
//    cpu_time += timer::stop(timer_dmi);
//    return -2. * state.material.D_atom/(constants::mu0 * state.material.p) * join(3, Hdx, Hdy, Hdz); // TODO FACTOR 2?
//    //print ("H", H);
//    //return H;
//}



////VERSION FOR fixed e_z=(0, 0, 1)
//#include "atomistic_dmi.hpp"

//using namespace af;
////Energy calculation
////E=-mu0/2 integral(M . H) dx
//float AtomisticDmiField::E(const State& state){
//  return - constants::mu0/2. * state.material.p * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3));
//  //return -constants::mu0/2. * state.material.p * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3));
//  //return -constants::mu0/2. * state.Ms * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
//}
//
//AtomisticDmiField::AtomisticDmiField (const Mesh& mesh, const Material& material){
//  //initialize finite difference first order derivative filter
//  float norm=sqrt(pow(material.D_axis[0], 2)+ pow(material.D_axis[1], 2) + pow(material.D_axis[2], 2));
//  n=array(mesh.n0, mesh.n1, mesh.n2, 3, f32);
//  n(span, span, span, 0)=material.D_axis[0]/norm;
//  n(span, span, span, 1)=material.D_axis[1]/norm;
//  n(span, span, span, 2)=material.D_axis[2]/norm;
//  //print("n", n);
//
//  //initialize finite difference first order derivative filter
//  filtr_atom_dmi_x=constant(0.0, 3, 3, 3, 3, f32);
//  filtr_atom_dmi_y=constant(0.0, 3, 3, 3, 3, f32);
//  filtr_atom_dmi_z=constant(0.0, 3, 3, 3, 3, f32);
//  //dmx/dx
//  filtr_atom_dmi_x(0, 1, 1, 2)= 1;
//  filtr_atom_dmi_x(2, 1, 1, 2)=-1;
//  //print("filtr", filtr_atom_dmi_x);
//
//  //dmy/dy
//  filtr_atom_dmi_y(1, 0, 1, 2)= 1;
//  filtr_atom_dmi_y(1, 2, 1, 2)=-1;
//  //print("filtr_y", filtr_atom_dmi_y);
//
//  //dmz/dz
//  filtr_atom_dmi_z(0, 1, 1, 0)=-1;
//  filtr_atom_dmi_z(2, 1, 1, 0)= 1;
//  filtr_atom_dmi_z(1, 0, 1, 1)=-1;
//  filtr_atom_dmi_z(1, 2, 1, 1)= 1;
//  //print("filtr_z", filtr_atom_dmi_z);
//}
//
//array AtomisticDmiField::h(const State& state){
//  timer_dmi = timer::start();
//
//  //print("m", state.m);
//  array convx = convolve(state.m, filtr_atom_dmi_x, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//  //print("convx", convx);
//  array convy = convolve(state.m, filtr_atom_dmi_y, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//  //print("convy", convy);
//  array convz = convolve(state.m, filtr_atom_dmi_z, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//  //print("convz", convz);
//  array res = join(3, convx(span, span, span, 2), convy(span, span, span, 2), convz(span, span, span, 0)+convz(span, span, span, 1));
//  //print("res", res);
//  if(state.material.afsync) sync();
//  cpu_time += timer::stop(timer_dmi);
//  return -state.material.D_atom/(2.*constants::mu0*state.material.p) * res;
//  //return -state.material.D/2. * res;
//}
}// namespace magnumaf
