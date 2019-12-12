#include "atomistic_exchange.hpp"
#include "../func.hpp"

namespace magnumafcpp{

using namespace af;

//Energy calculation
//Eex=-mu0/2 integral(M . Hex) dx

double AtomisticExchangeField::E(const State& state){
  return -constants::mu0/2. *state.Ms * afvalue(sum(sum(sum(sum(h(state)*state.m, 0), 1), 2), 3));
}

double AtomisticExchangeField::E(const State& state, const af::array& h){
  return -constants::mu0/2. *state.Ms * afvalue(sum(sum(sum(sum(h * state.m, 0), 1), 2), 3));
}

AtomisticExchangeField::AtomisticExchangeField (double J_atom) : J_atom(J_atom)
{
}

array AtomisticExchangeField::h(const State& state){

    array filtr=constant(0.0, 3, 3, 3, 3, f64);
    filtr(0, 1, 1, span)= 1.;
    filtr(2, 1, 1, span)= 1.;
    filtr(1, 0, 1, span)= 1.;
    filtr(1, 2, 1, span)= 1.;
    filtr(1, 1, 0, span)= 1.;
    filtr(1, 1, 2, span)= 1.;
    //if(state.material.hexagonal_close_packed == true){
    //    std::cout << "WARNING: Experimental hcp exchange" << std::endl;
    //    filtr(0, 1, 1, span)= 1.;
    //    filtr(0, 2, 1, span)= 1.;// hex lattice
    //    filtr(2, 1, 1, span)= 1.;
    //    filtr(2, 0, 1, span)= 1.;// hex lattice
    //    filtr(1, 0, 1, span)= 1.;
    //    filtr(1, 2, 1, span)= 1.;
    //    filtr(1, 1, 0, span)= 1.;
    //    filtr(1, 1, 2, span)= 1.;
    //    //TODO hex in z: filtr(, , , span)= 1.;// hex lattice
    //    //NOTE: numers of NN is 12, so 3 in +z and 3 in +z slice
    //    af::print("filtr", filtr);

    //}
    //else if(state.material.atom_fcc=true){
    //    //https://www.physics-in-a-nutshell.com/article/11/close-packed-structures-fcc-and-hcp
    //}
    //else{
    //}

  timer_solve = timer::start();
  //convolution
  array mj = convolve(state.m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);

  if(state.afsync) af::sync();
  cpu_time += timer::stop(timer_solve);
  return J_atom/(constants::mu0*state.Ms)* mj;
}
  //return state.material.J/(state.mesh.dx*constants::mu0) * mj;
  //return state.material.J/(state.mesh.dx*constants::mu0*state.Ms) * mj;
  //return state.material.J/state.mesh.dx * mj;




//  //testing
//  array m = constant(0.0, state.mesh.dims, f64);
//  m(span, span, span, 2)=1;
//  //m(span, span, span, 1)=2;
//  //m(span, span, span, 2)=3;
//  print("m", m);
//  array filtr_x=constant(0.0, 3, 3, 3, 3, f64);
//  print("test", convolve(m, filtr_x, AF_CONV_DEFAULT, AF_CONV_SPATIAL));









  //return -state.material.J/(2. * state.mesh.dx) * mj;//TODO not dx
  //array result =-state.material.J/(2. * state.mesh.dx) * mj;
  //Wrong: this would be the Energy, not the field:   array result =tile(-state.material.J/2. * sum(state.m*mj, 3), 1, 1, 1, 3);
  //return result;
  //return  -state.material.J/2. * sum(state.m*mj, 3);



//AtomisticExchangeField::AtomisticExchangeField (Mesh meshin, Material paramin) : material(paramin), mesh(meshin){
//  array exch = array(mesh.n0, mesh.n1, mesh.n2, 3, f64);
//  //initialize filters
//  filtr=constant(0.0, 3, 3, 3, f64);
//  filtr(1, 1, 1)= -6 / (pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//
//  filtr(0, 1, 1)= 1 / pow(mesh.dx, 2);
//  filtr(2, 1, 1)= 1 / pow(mesh.dx, 2);
//
//  filtr(1, 0, 1)= 1 / pow(mesh.dy, 2);
//  filtr(1, 2, 1)= 1 / pow(mesh.dy, 2);
//
//  filtr(1, 1, 0)= 1 / pow(mesh.dz, 2);
//  filtr(1, 1, 2)= 1 / pow(mesh.dz, 2);
//}
//array AtomisticExchangeField::solve(array m){
//  timer_exchsolve = timer::start();
//
//  timer_conv = timer::start();
//
//  //convolution
//  //exch = convolve(m, filtr);
//  exch = convolve(m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//
//  if(state.material.afsync) sync();
//  time_conv += timer::stop(timer_conv);
//
//  //Accounting for boundary conditions by adding initial m values on the boundaries by adding all 6 boundary surfaces
//  timer_edges = timer::start();
//  exch(0, span, span, span)+=m(0 , span, span, span)/ pow(mesh.dx, 2);
//  exch(-1, span, span, span)+=m(-1, span, span, span)/ pow(mesh.dx, 2);
//
//
//  exch(span, 0 , span, span)+=m(span, 0 , span, span)/ pow(mesh.dy, 2);
//  exch(span, -1, span, span)+=m(span, -1, span, span)/ pow(mesh.dy, 2);
//
//  exch(span, span, 0 , span)+=m(span, span, 0 , span)/ pow(mesh.dz, 2);
//  exch(span, span, -1, span)+=m(span, span, -1, span)/ pow(mesh.dz, 2);
//  if(state.material.afsync) sync();
//  time_edges += timer::stop(timer_edges);
//  time_exchsolve += timer::stop(timer_exchsolve);
//  return  (2.* material.A)/(constants::mu0*state.Ms) * exch;
//}
//
//void showdims2(const array& a){
//  std::cout<<"Exchange matrix: dims="<<a.dims(0)<<"\t"<<a.dims(1)<<"\t"<<a.dims(2)<<"\t"<<a.dims(3)<<std::endl;
//}
//
////Function returns index
//int AtomisticExchangeField::findex(int i0, int i1, int i2, int im, int id){
//  return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
//}
//
////Inner index (index per matrix column)
//int AtomisticExchangeField::findex(int i0, int i1, int i2, int im){
//  return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*im));
//}
//
//AtomisticExchangeField::AtomisticExchangeField (Mesh meshin, Material paramin) : material(paramin), mesh(meshin){
//  const int dimension=mesh.n0*mesh.n1*mesh.n2*3;
//  double* vmatr = NULL;
//  vmatr = new double[dimension*dimension];
//  for (int id = 0; id < dimension; id++){
//    for (int im = 0; im < 3; im++){
//      for (int i2 = 0; i2 < mesh.n2; i2++){
//        for (int i1 = 0; i1 < mesh.n1; i1++){
//          for (int i0 = 0; i0 < mesh.n0; i0++){
//            const int index=i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
//            vmatr[index]=0.;
//          }
//        }
//      }
//    }
//  }
//
//  for (int id = 0; id < dimension; id++){
//    for (int im = 0; im < 3; im++){
//      for (int i2 = 0; i2 < mesh.n2; i2++){
//        for (int i1 = 0; i1 < mesh.n1; i1++){
//          for (int i0 = 0; i0 < mesh.n0; i0++){
//            const int ind=findex(i0, i1, i2, im);
//            //const int index=i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
//            //const int ind=i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*im));
//            if(ind==id) {
//              vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//              //x
//              if(i0==0){
//                vmatr[findex(i0  , i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                if (mesh.n0>1) vmatr[findex(i0+1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//              }
//              if (i0==mesh.n0-1){
//                vmatr[findex(i0  , i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                if (mesh.n0>1) vmatr[findex(i0-1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//              }
//              if(i0>0 && i0< mesh.n0-1){
//                vmatr[findex(i0-1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                vmatr[findex(i0+1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//              }
//
//              //y
//              if(i1==0){
//                vmatr[findex(i0, i1  , i2, im, id)]+= 1./pow(mesh.dy, 2);
//                if (mesh.n1>1) vmatr[findex(i0, i1+1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//              }
//              if (i1==mesh.n1-1){
//                vmatr[findex(i0, i1  , i2, im, id)]+= 1./pow(mesh.dy, 2);
//                if (mesh.n1>1) vmatr[findex(i0, i1-1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//              }
//              if(i1>0 && i1< mesh.n1-1){
//                vmatr[findex(i0, i1-1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//                vmatr[findex(i0, i1+1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//              }
//
//              //z
//              if (i2==0){
//                vmatr[findex(i0, i1, i2  , im, id)]+= 1./pow(mesh.dz, 2);
//                if (mesh.n2>1) vmatr[findex(i0, i1, i2+1, im, id)]+= 1./pow(mesh.dz, 2);
//              }
//              if (i2==mesh.n2-1){
//                vmatr[findex(i0, i1, i2  , im, id)]+= 1./pow(mesh.dz, 2);
//                if (mesh.n2>1) vmatr[findex(i0, i1, i2-1, im, id)]+= 1./pow(mesh.dz, 2);
//              }
//              if(i2>0 && i2< mesh.n2-1){
//                vmatr[findex(i0, i1, i2-1, im, id)]+= 1./pow(mesh.dz, 2);
//                vmatr[findex(i0, i1, i2+1, im, id)]+= 1./pow(mesh.dz, 2);
//              }
//
//            //  if(i1==0){
//            //  }
//            //  else if (i1==mesh.n1-1){
//            //  }
//            //  else
//            //    vmatr[findex(i0, i1-1, i2, im, id)]+= 1.;
//            //    vmatr[findex(i0, i1+1, i2, im, id)]+= 1.;
//            //  if(i2==0){
//            //  }
//            //  else if (i2==mesh.n2-1){
//            //  }
//            //  else
//            //    vmatr[findex(i0, i1, i2-1, im, id)]+= 1.;
//            //    vmatr[findex(i0, i1, i2+1, im, id)]+= 1.;
//
//              //if((index > 0) && (index < dimension*dimension -1 )){
//              //  vmatr[findex(i0, i1-1, i2, im, id)]+= 1.;
//              //  vmatr[findex(i0, i1+1, i2, im, id)]+= 1.;
//              //}
//            }
//            //if(i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im)))==id) vmatr[index]=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//          }
//        }
//      }
//    }
//  }
//  array fullmatr(dimension, dimension, vmatr);
//  delete [] vmatr;
//  vmatr = NULL;
//  showdims2(fullmatr);
//  //print("fullmatr", fullmatr);
////  print("fullmatr", fullmatr(span, 0, 0, 0));
////  print("fullmatr", fullmatr(span, 1, 0, 0));
////  print("fullmatr", fullmatr(span, -4, 0, 0));
////  print("fullmatr", fullmatr(span, 30, 0, 0));
////  print("fullmatr", fullmatr(span, -3, 0, 0));
////  print("fullmatr", fullmatr(span, -1, 0, 0));
//  matr=sparse(fullmatr);
//
//  std::cout << "Sparsity of matr = "
//            << (float)sparseGetNNZ(matr) / (float)matr.elements()
//            << std::endl;
//}
//
//
//
//
//array AtomisticExchangeField::solve(array m){
////print("m", flat(m));
//  timer_exchsolve = timer::start();
//  array exch = matmul(matr, flat(m));
////print("ex0", exch);
//  exch=moddims(exch, mesh.n0, mesh.n1, mesh.n2, 3);
////print("ex1", flat(exch));
//
//  exch.eval();
//  if(state.material.afsync) sync();
//  time_exchsolve += timer::stop(timer_exchsolve);
//
//  return  (2.* material.A)/(constants::mu0*state.Ms) * exch;
//}
//
//
}// namespace magnumafcpp
