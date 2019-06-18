#include "micro_dmi.hpp"
using namespace af;
void showdims(const array& a){
  std::cout<<"Exchange matrix: dims="<<a.dims(0)<<"\t"<<a.dims(1)<<"\t"<<a.dims(2)<<"\t"<<a.dims(3)<<std::endl;
}
void apply_boundary_condition(array& hfield, const State& state);
//Energy calculation
//E=-mu0/2 integral(M . H) dx
double DmiField::E(const State& state){
  return -constants::mu0/2. * state.Ms * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3)) * mesh.dx * mesh.dy * mesh.dz;
}

double DmiField::E(const State& state, const af::array& h){
  return -constants::mu0/2. * state.Ms * afvalue(sum(sum(sum(sum(h * state.m, 0), 1), 2), 3)) * mesh.dx * mesh.dy * mesh.dz;
}

DmiField::DmiField (Mesh meshin, Material paramin) : material(paramin), mesh(meshin){
  //Normal vector
  double norm=sqrt(pow(material.D_axis[0], 2)+ pow(material.D_axis[1], 2) + pow(material.D_axis[2], 2));
  n=array(mesh.n0, mesh.n1, mesh.n2, 3, f64);
  n(span, span, span, 0)=material.D_axis[0]/norm;
  n(span, span, span, 1)=material.D_axis[1]/norm;
  n(span, span, span, 2)=material.D_axis[2]/norm;
  //print("n", n);

  //initialize finite difference first order derivative filter
  //Central finite difference in 2nd order
  //https://en.wikipedia.org/wiki/Finite_difference_coefficient
  filtr_fd1=constant(0.0, 3, 3, 3, 3, f64);
  //dmx/dx
  filtr_fd1(0, 1, 1, 0)= 1 / (2.*mesh.dx);
  filtr_fd1(2, 1, 1, 0)=-1 / (2.*mesh.dx);

  //dmy/dy
  filtr_fd1(1, 0, 1, 1)= 1 / (2.*mesh.dy);
  filtr_fd1(1, 2, 1, 1)=-1 / (2.*mesh.dy);

  //dmz/dz
  filtr_fd1(1, 1, 0, 2)= 1 / (2.*mesh.dz);
  filtr_fd1(1, 1, 2, 2)=-1 / (2.*mesh.dz);
}


array DmiField::h(const State& state){
  timer_dmi = timer::start();

  //First: n(div m)
  //Gradient and edges
  array first = convolve(state.m, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);

  //correct_edges(first, state.m);
  apply_boundary_condition(first, state);
  //Make Divergence
  first=sum(first, 3);
  first=tile(first, 1, 1, 1, 3);
  first=n*first;

  //Second: fd1(n . m)
  //Dot product
  array second = sum(n*state.m, 3);
  //Expand for fd1 convolution
  second = tile(second, 1, 1, 1, 3);
  second = convolve(second, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
  //correct_edges(second, state.m);
  apply_boundary_condition(second, state);

  if(state.material.afsync) af::sync();
  cpu_time += timer::stop(timer_dmi);
  return 2.* material.D/(constants::mu0*state.Ms) * (first-second);//Note: Js=mu0*Ms
}


//EXPERIMENTAL
//This version assumes the same boundaries as in the exchange filed of 70Lines of Numpy
//On the edges, we assume (not existing) outer cells with the same values as the boundary cells
// in the middle ...|-1/2|  0 |1/2|...
// on the edges         ||  0 |1/2|...
// we want         |-1/2||  0 |1/2|...
// thus we take         ||-1/2|1/2|...
// so after the convolution we have to add the edges with -1/(2*dx)

void DmiField::correct_edges(array& out, const array& in){
  //Lower x edge:
  out( 0, span, span, 0)+= -0.5* in( 0, span, span, 0)/mesh.dx;
  //Upper x edge:
  out(-1, span, span, 0)-= -0.5* in(-1, span, span, 0)/mesh.dx;

  //Lower y edge:
  out(span, 0, span, 1)+= -0.5* in(span, 0, span, 1)/mesh.dy;
  //Upper y edge:
  out(span, -1, span, 1)-= -0.5* in(span, -1, span, 1)/mesh.dy;

  //z
  if(in.dims(2)==1){
    out(span, span, span, 2)=0.;
  }
  else{
    //Lower z edge:
    out(span, span, 0, 2)+= -0.5* in(span, span, 0, 2)/mesh.dz;
    //Upper z edge:
    out(span, span, -1, 2)-= -0.5* in(span, span, -1, 2)/mesh.dz;
  }
}

//--------------------------------------------------------------------------------------------------------------------------------
// EXPERIMENTAL
// TODO unit tests !
// Boundary Conditions according to DOI: 10.1103/PhysRevB.88.184422 Skyrmion confinement in ultrathin film nanostructures ...
// dm/dn = 1/xi (n_DM x n_surface) x m ; xi = 2 A/D
// here:
// dm/dn = 1/xi (D_axis x n_surface) x m ; xi = 2 A/D
void apply_boundary_condition(array& hfield, const State& state){
    //DM Vector:
    const array n_DM(1, 1, 1, 3, state.material.D_axis);
    double A = 0;//TODO set exchange A as class member and pass in constructor

    if(state.m.dims(0)==1){
        hfield(span, span, span, 0)=0.;
    }
    else{
        // low x boundary:
        // n_surface=(-1, 0, 0)
        // dm/dn=dm/d(-x)= (m_-1 - m_1) / 2 * dx = 1/xi (D_axis x n_surface) x m
        // => m_-1 = m_1 + 2 * dx * 1/xi (D_axis x n_surface) x m_0
        const double c_n_x_surface_low[]={-1, 0, 0};
        const array n_x_surface_low(1, 1, 1, 3, c_n_x_surface_low); // normal vector to the x-surface at boundary with lower index i=0
        const array n_DMxn_x_surf_low=tile(cross4(n_DM, n_x_surface_low), 1, state.m.dims(1), state.m.dims(2), 1);
        const array x_minus_1 = state.m(1, span, span, 0) + 2 * state.mesh.dx * (state.material.D/(2*A))  * cross4(n_DMxn_x_surf_low , state.m(0, span, span, span))(span, span, span, 0);
        hfield(0, span, span, 0)+= - 0.5 * x_minus_1 / state.mesh.dx; // Minus due to: (m_{i+1} - m_{i-1})/( 2*dx )  with m_{i-1} being replaced

        // high x boundary:
        // n_surface=(1, 0, 0)
        // dm/dn=dm/d(x)= (m_{n+1} - m_{n-1}) / 2 * dx = 1/xi (D_axis x n_surface) x m_n
        // => m_{i+1} = m_{i-1} + 2 * dx * 1/xi (D_axis x n_surface) x m_i
        const double c_n_x_surface_high[]={1, 0, 0};
        const array n_x_surface_high(1, 1, 1, 3, c_n_x_surface_high); // normal vector to the x-surface at boundary with higher index i=0
        const array n_DMxn_x_surf_high=tile(cross4(n_DM, n_x_surface_high), 1, state.m.dims(1), state.m.dims(2), 1);
        const array x_i_plus_1 = state.m(-1, span, span, 0) + 2 * state.mesh.dx * (state.material.D/(2*A))  * cross4(n_DMxn_x_surf_high , state.m(-1, span, span, span))(span, span, span, 0);
        hfield(-1, span, span, 0)+= 0.5 * x_i_plus_1 / state.mesh.dx;
    }

    if(state.m.dims(1)==1){
        hfield(span, span, span, 1)=0.;
    }
    else{
        // low y boundary:
        // n_surface=(0, -1, 0)
        // dm/dn=dm/d(-y)= (m_-1 - m_1) / 2 * dx = 1/xi (D_axis x n_surface) x m_0
        // => m_-1 = m_1 + 2 * dx * 1/xi (D_axis x n_surface) x m_0
        const double c_n_y_surface_low[]={0, -1, 0};
        const array n_y_surface_low(1, 1, 1, 3, c_n_y_surface_low); // normal vector to the y-surface at boundary with lower indey i=0
        const array n_DMxn_y_surf_low=tile(cross4(n_DM, n_y_surface_low), state.m.dims(0), 1, state.m.dims(2), 1);
        const array y_minus_1 = state.m(span, 1, span, 1) + 2 * state.mesh.dy * (state.material.D/(2*A))  * cross4(n_DMxn_y_surf_low , state.m(span, 0, span, span))(span, span, span, 1);
        hfield(span, 0, span, 1)+= - 0.5 * y_minus_1 / state.mesh.dy; // Minus due to: (m_{i+1} - m_{i-1})/( 2*dy )  with m_{i-1} being replaced

        // high y boundary:
        // n_surface=(0, 1, 0)
        const double c_n_y_surface_high[]={0, 1, 0};
        const array n_y_surface_high(1, 1, 1, 3, c_n_y_surface_high); // normal vector to the y-surface at boundary with higher indey i=0
        const array n_DMxn_y_surf_high=tile(cross4(n_DM, n_y_surface_high), state.m.dims(0), 1, state.m.dims(2), 1);
        const array y_i_plus_1 = state.m(span, -1, span, 1) + 2 * state.mesh.dy * (state.material.D/(2*A))  * cross4(n_DMxn_y_surf_high , state.m(span, -1, span, span))(span, span, span, 1);
        hfield(span, -1, span, 1)+= 0.5 * y_i_plus_1 / state.mesh.dy; // Minus due to: (m_{i+1} - m_{i-1})/( 2*dy )  with m_{i-1} being replaced
    }

    //z
    if(state.m.dims(2)==1){
      hfield(span, span, span, 2)=0.;
    }
    else{
        // low z boundary:
        // n_surface=(0, 0, -1)
        // dm/dn=dm/d(-z)= (m_-1 - m_1) / 2 * dz = 1/xi (D_axis x n_surface) x m
        // => m_-1 = m_1 + 2 * dz * 1/xi (D_axis x n_surface) x m_0
        const double c_n_z_surface_low[]={0, 0, -1};
        const array n_z_surface_low(1, 1, 1, 3, c_n_z_surface_low); // normal vector to the z-surface at boundary with lower index i=0
        const array n_DMxn_z_surf_low=tile(cross4(n_DM, n_z_surface_low), state.m.dims(0), state.m.dims(1), 1, 1);
        const array z_minus_1 = state.m(span, span, 1, 2) + 2 * state.mesh.dz * (state.material.D/(2*A))  * cross4(n_DMxn_z_surf_low , state.m(span, span, 0, span))(span, span, span, 2);
        hfield(span, span, 0, 2)+= - 0.5 * z_minus_1 / state.mesh.dz; // Minus due to: (m_{i+1} - m_{i-1})/( 2*dz )  with m_{i-1} being replaced

        // high z boundary:
        // n_surface=(1, 0, 0)
        // dm/dn=dm/d(z)= (m_{n+1} - m_{n-1}) / 2 * dz = 1/xi (D_axis x n_surface) x m_n
        // => m_{i+1} = m_{i-1} + 2 * dz * 1/xi (D_axis x n_surface) x m_i
        const double c_n_z_surface_high[]={0, 0, 1};
        const array n_z_surface_high(1, 1, 1, 3, c_n_z_surface_high); // normal vector to the z-surface at boundary with higher index i=0
        const array n_DMxn_z_surf_high=tile(cross4(n_DM, n_z_surface_high), state.m.dims(0), state.m.dims(1), 1, 1);
        const array z_i_plus_1 = state.m(span, span, -1, 2) + 2 * state.mesh.dz * (state.material.D/(2*A))  * cross4(n_DMxn_z_surf_high , state.m(span, span, -1, span))(span, span, span, 2);
        hfield(span, span, -1, 2)+= 0.5 * z_i_plus_1 / state.mesh.dz;
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

//void DmiField::correct_edges(array& out, const array& in){
//  //Lower x bound: after convolve it is:  1/2 * (1)
//  //i.e.                  filtr_grad(0, 1, 1, 0) * in(0, span, span, 0)
//  //We want it to be:                  -1*(0) + 1* (1)
//  //So we take                        1/2*(1) - 1* (0)
//
//  //Lower x edge:
//  out( 0, span, span, 0)+= in( 1, span, span, 0)/(2.*mesh.dx) - in( 0, span, span, 0)/mesh.dx;
//  //Upper x edge:
//  out(-1, span, span, 0)-= in(-2, span, span, 0)/(2.*mesh.dx) - in(-1, span, span, 0)/mesh.dx;
//
//  //Lower y edge:
//  out(span, 0, span, 1)+= in(span, 1, span, 1)/(2.*mesh.dy) - in(span, 0, span, 1)/mesh.dy;
//  //Upper y edge:
//  out(span, -1, span, 1)-= in(span, -2, span, 1)/(2.*mesh.dy) - in(span, -1, span, 1)/mesh.dy;
//
//  //z
//  if(in.dims(2)==1){
//    out(span, span, span, 2)=0.;
//  }
//  else{
//    //Lower z edge:
//    out(span, span, 0, 2)+= in(span, span, 1, 2)/(2.*mesh.dz) - in(span, span, 0, 2)/mesh.dz;
//    //Upper z edge:
//    out(span, span, -1, 2)-= in(span, span, -2, 2)/(2.*mesh.dz) - in(span, span, -1, 2)/mesh.dz;
//  }
//}















//DmiField::DmiField (Mesh meshin, Material paramin) : material(paramin), mesh(meshin){
//  //Normal vector
//  n=array(mesh.n0, mesh.n1, mesh.n2, 3, f64);
//  n(span, span, span, 0)=nx;
//  n(span, span, span, 1)=ny;
//  n(span, span, span, 2)=nz;
//  //print("n", n);
//
//  //TODO:
//  //initialize filters
//  filtr=constant(0.0, 3, 3, 3, f64);
//  //filtr(1, 1, 1)= -6 / (pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//
//  filtr(0, 1, 1)=-1 / (2.*mesh.dx);
//  filtr(2, 1, 1)= 1 / (2.*mesh.dx);
//
//  filtr(1, 0, 1)=-1 / (2.*mesh.dy);
//  filtr(1, 2, 1)= 1 / (2.*mesh.dy);
//
//  filtr(1, 1, 0)=-1 / (2.*mesh.dz);
//  filtr(1, 1, 2)= 1 / (2.*mesh.dz);
//}
//array DmiField::h(array m){
//  timer_dmi = timer::start();
//  array first = convolve(m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//  showdims(first);
//  first=n*tile(sum(first, 3), 1, 1, 1, 3);
//  showdims(first);
//  array second = convolve(tile(sum(n*m, 3), 1, 1, 1, 3), filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//  showdims(second);
//
//  if(state.material.afsync) sync();
//  cpu_time += timer::stop(timer_dmi);
//  return  -2.* material.D/material.Js * (first-second);//TODO Js not set
//}
