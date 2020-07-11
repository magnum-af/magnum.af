#include "micro_exch.hpp"
#include "../func.hpp"

namespace magnumafcpp {

// Energy calculation
// Eex=-mu0/2 integral(M . Hex) dx
// Virtual method is overwritten as to use h_withedges
// Note: maybe this is irrelevant and can be dropped.
double ExchangeField::E(const State& state) {
    if (state.Ms_field.isempty()) {
        return -constants::mu0 / 2. * state.Ms *
               af::sum(af::sum(af::sum(af::sum(h_withedges(state) * state.m, 0), 1), 2), 3).scalar<double>() *
               state.mesh.dx * state.mesh.dy * state.mesh.dz;
    } else {
        return -constants::mu0 / 2. *
               af::sum(af::sum(af::sum(af::sum(state.Ms_field * h_withedges(state) * state.m, 0), 1), 2), 3)
                   .scalar<double>() *
               state.mesh.dx * state.mesh.dy * state.mesh.dz;
    }
}

/// Constructor for a global exchange constant
ExchangeField::ExchangeField(double A) : A(A) {}

/// Constructor taking spacially varying exchange constant af af::array
ExchangeField::ExchangeField(af::array A_field)
    : A_field(A_field.dims(3) == 1 ? af::tile(A_field, 1, 1, 1, 3) : A_field) {
    printf("\33[1;31mWarning:\33[0m ExchangeField: This is legacy code, to use "
           "spacially varying A values with correct jump conditions, use "
           "SparseExchangeField or RKKYExchangeField!\n");
}

// For wrapping only
ExchangeField::ExchangeField(long int A_field_ptr)
    : A_field((*(new af::array(*((void**)A_field_ptr)))).dims(3) == 1
                  ? af::tile(*(new af::array(*((void**)A_field_ptr))), 1, 1, 1, 3)
                  : *(new af::array(*((void**)A_field_ptr)))) {
    printf("\33[1;31mWarning:\33[0m ExchangeField: This is legacy code, to use "
           "spacially varying A values with correct jump conditions, use "
           "SparseExchangeField or RKKYExchangeField!\n");
}

af::array ExchangeField::h_withedges(const State& state) {
    timer_exchsolve = af::timer::start();
    af::array filtr = af::constant(0.0, 3, 3, 3, f64);
    // Note: skipped as this term falls out int cross product: //filtr(1, 1, 1)=
    // -6 / (pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
    filtr(0, 1, 1) = 1 / pow(state.mesh.dx, 2);
    filtr(2, 1, 1) = 1 / pow(state.mesh.dx, 2);
    filtr(1, 0, 1) = 1 / pow(state.mesh.dy, 2);
    filtr(1, 2, 1) = 1 / pow(state.mesh.dy, 2);
    filtr(1, 1, 0) = 1 / pow(state.mesh.dz, 2);
    filtr(1, 1, 2) = 1 / pow(state.mesh.dz, 2);
    // Convolution
    af::array exch = convolve(state.m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);

    // Accounting for boundary conditions by adding initial m values on the
    // boundaries by adding all 6 boundary surfaces
    timer_edges = af::timer::start();
    exch(0, af::span, af::span, af::span) += state.m(0, af::span, af::span, af::span) / pow(state.mesh.dx, 2);
    exch(-1, af::span, af::span, af::span) += state.m(-1, af::span, af::span, af::span) / pow(state.mesh.dx, 2);

    exch(af::span, 0, af::span, af::span) += state.m(af::span, 0, af::span, af::span) / pow(state.mesh.dy, 2);
    exch(af::span, -1, af::span, af::span) += state.m(af::span, -1, af::span, af::span) / pow(state.mesh.dy, 2);

    exch(af::span, af::span, 0, af::span) += state.m(af::span, af::span, 0, af::span) / pow(state.mesh.dz, 2);
    exch(af::span, af::span, -1, af::span) += state.m(af::span, af::span, -1, af::span) / pow(state.mesh.dz, 2);

    if (state.afsync)
        af::sync();
    time_edges += af::timer::stop(timer_edges);
    computation_time_heff += af::timer::stop(timer_exchsolve);
    if (state.Ms_field.isempty() && this->A_field.isempty()) {
        return (2. * this->A) / (constants::mu0 * state.Ms) * exch;
    } else if (!state.Ms_field.isempty() && this->A_field.isempty()) {
        af::array heff = (2. * this->A) / (constants::mu0 * state.Ms_field) * exch;
        replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    } else if (state.Ms_field.isempty() && !this->A_field.isempty()) {
        return (2. * this->A_field) / (constants::mu0 * state.Ms) * exch;
    } else {
        af::array heff = (2. * this->A_field) / (constants::mu0 * state.Ms_field) * exch;
        replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// Terms proportional to m dorp out in the cross product of the LLG and thus is
// neglected as arrayfire is extremely slow with indexing operations NOTE: This
// yields no longer the physical exchange field but optimizes the caluclation
af::array ExchangeField::h(const State& state) {
    timer_exchsolve = af::timer::start();
    af::array filtr = af::constant(0.0, 3, 3, 3, f64);
    // Note: skipped as this term falls out int cross product: //filtr(1, 1, 1)=
    // -6 / (pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
    filtr(0, 1, 1) = 1 / pow(state.mesh.dx, 2);
    filtr(2, 1, 1) = 1 / pow(state.mesh.dx, 2);
    filtr(1, 0, 1) = 1 / pow(state.mesh.dy, 2);
    filtr(1, 2, 1) = 1 / pow(state.mesh.dy, 2);
    filtr(1, 1, 0) = 1 / pow(state.mesh.dz, 2);
    filtr(1, 1, 2) = 1 / pow(state.mesh.dz, 2);
    af::array exch = convolve(state.m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    if (state.afsync)
        af::sync();
    computation_time_heff += af::timer::stop(timer_exchsolve);
    if (state.Ms_field.isempty() && this->A_field.isempty()) {
        return (2. * this->A) / (constants::mu0 * state.Ms) * exch;
    } else if (!state.Ms_field.isempty() && this->A_field.isempty()) {
        af::array heff = (2. * this->A) / (constants::mu0 * state.Ms_field) * exch;
        replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    } else if (state.Ms_field.isempty() && !this->A_field.isempty()) {
        return (2. * this->A_field) / (constants::mu0 * state.Ms) * exch;
    } else {
        af::array heff = (2. * this->A_field) / (constants::mu0 * state.Ms_field) * exch;
        replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// void showdims2(const array& a){
//  std::cout<<"Exchange matrix:
//  dims="<<a.dims(0)<<"\t"<<a.dims(1)<<"\t"<<a.dims(2)<<"\t"<<a.dims(3)<<std::endl;
//}

////Version with switch conv/sparseMatMul dropping the edges
//#include "exch.hpp"
//
// void showdims2(const array& a){
//  std::cout<<"Exchange matrix:
//  dims="<<a.dims(0)<<"\t"<<a.dims(1)<<"\t"<<a.dims(2)<<"\t"<<a.dims(3)<<std::endl;
//}
//
////Energy calculation
////Eex=-mu0/2 integral(M . Hex) dx
// double ExchangeField::E(const State& state){
//  return -constants::mu0/2. * state.Ms *
//  afvalue(af::sum(af::sum(af::sum(af::sum(h(state)*state.m, 0), 1), 2), 3)) *
//  mesh.dx * mesh.dy * mesh.dz;
//}
//
////Function returns index
// int ExchangeField::findex(int i0, int i1, int i2, int im, int id){
//  return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
//}
//
////Inner index (index per matrix column)
// int ExchangeField::findex(int i0, int i1, int i2, int im){
//  return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*im));
//}
//
// ExchangeField::ExchangeField (Mesh meshin, Material paramin) :
// material(paramin), mesh(meshin){
//  if(mesh.n0*mesh.n1*mesh.n2>8128){
//    //initialize filters
//    filtr=constant(0.0, 3, 3, 3, f64);
//    //filtr(1, 1, 1)= -6 / (pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//
//    filtr(0, 1, 1)= 1 / pow(mesh.dx, 2);
//    filtr(2, 1, 1)= 1 / pow(mesh.dx, 2);
//
//    filtr(1, 0, 1)= 1 / pow(mesh.dy, 2);
//    filtr(1, 2, 1)= 1 / pow(mesh.dy, 2);
//
//    filtr(1, 1, 0)= 1 / pow(mesh.dz, 2);
//    filtr(1, 1, 2)= 1 / pow(mesh.dz, 2);
//  }
//  //Currently better performance for small systems with matmul
//  else{
//    const int dimension=mesh.n0*mesh.n1*mesh.n2*3;
//    double* vmatr = NULL;
//    vmatr = new double[dimension*dimension];
//    for (int id = 0; id < dimension; id++){
//      for (int im = 0; im < 3; im++){
//        for (int i2 = 0; i2 < mesh.n2; i2++){
//          for (int i1 = 0; i1 < mesh.n1; i1++){
//            for (int i0 = 0; i0 < mesh.n0; i0++){
//              const int index=i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
//              vmatr[index]=0.;
//            }
//          }
//        }
//      }
//    }
//
//    for (int id = 0; id < dimension; id++){
//      for (int im = 0; im < 3; im++){
//        for (int i2 = 0; i2 < mesh.n2; i2++){
//          for (int i1 = 0; i1 < mesh.n1; i1++){
//            for (int i0 = 0; i0 < mesh.n0; i0++){
//              const int ind=findex(i0, i1, i2, im);
//              if(ind==id) {
//                //vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx,
//                2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//                //x
//                if(i0==0){
//                  //vmatr[findex(i0  , i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                  if (mesh.n0>1) vmatr[findex(i0+1, i1, i2, im,
//                  id)]+= 1./pow(mesh.dx, 2);
//                }
//                if (i0==mesh.n0-1){
//                  //vmatr[findex(i0  , i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                  if (mesh.n0>1) vmatr[findex(i0-1, i1, i2, im,
//                  id)]+= 1./pow(mesh.dx, 2);
//                }
//                if(i0>0 && i0< mesh.n0-1){
//                  vmatr[findex(i0-1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                  vmatr[findex(i0+1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                }
//
//                //y
//                if(i1==0){
//                  //vmatr[findex(i0, i1  , i2, im, id)]+= 1./pow(mesh.dy, 2);
//                  if (mesh.n1>1) vmatr[findex(i0, i1+1, i2, im,
//                  id)]+= 1./pow(mesh.dy, 2);
//                }
//                if (i1==mesh.n1-1){
//                  //vmatr[findex(i0, i1  , i2, im, id)]+= 1./pow(mesh.dy, 2);
//                  if (mesh.n1>1) vmatr[findex(i0, i1-1, i2, im,
//                  id)]+= 1./pow(mesh.dy, 2);
//                }
//                if(i1>0 && i1< mesh.n1-1){
//                  vmatr[findex(i0, i1-1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//                  vmatr[findex(i0, i1+1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//                }
//
//                //z
//                if (i2==0){
//                  //vmatr[findex(i0, i1, i2  , im, id)]+= 1./pow(mesh.dz, 2);
//                  if (mesh.n2>1) vmatr[findex(i0, i1, i2+1, im,
//                  id)]+= 1./pow(mesh.dz, 2);
//                }
//                if (i2==mesh.n2-1){
//                  //vmatr[findex(i0, i1, i2  , im, id)]+= 1./pow(mesh.dz, 2);
//                  if (mesh.n2>1) vmatr[findex(i0, i1, i2-1, im,
//                  id)]+= 1./pow(mesh.dz, 2);
//                }
//                if(i2>0 && i2< mesh.n2-1){
//                  vmatr[findex(i0, i1, i2-1, im, id)]+= 1./pow(mesh.dz, 2);
//                  vmatr[findex(i0, i1, i2+1, im, id)]+= 1./pow(mesh.dz, 2);
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//    array fullmatr(dimension, dimension, vmatr);
//    delete [] vmatr;
//    vmatr = NULL;
//    //showdims2(fullmatr);
//    matr=sparse(fullmatr);
//
//    std::cout << "Sparsity of matr = "
//              << (float)sparseGetNNZ(matr) / (float)matr.elements()
//              << std::endl;
//    }
//}
//
// array ExchangeField::h(const State& state){
//  timer_exchsolve = af::timer::start();
//
//  if(mesh.n0*mesh.n1*mesh.n2>8128){
//    timer_conv = af::timer::start();
//    //convolution
//    array exch = convolve(state.m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//
//    if(state.material.afsync) sync();
//    time_conv += af::timer::stop(timer_conv);
//
//    //Accounting for boundary conditions by adding initial m values on the
//    boundaries by adding all 6 boundary surfaces
//    //timer_edges = af::timer::start();
//    //exch(0, af::span, af::span, af::span)+=state.m(0 , af::span, af::span,
//    af::span)/ pow(mesh.dx, 2);
//    //exch(-1, af::span, af::span, af::span)+=state.m(-1, af::span, af::span,
//    af::span)/ pow(mesh.dx, 2);
//    //
//    //
//    //exch(af::span, 0 , af::span, af::span)+=state.m(af::span, 0 , af::span,
//    af::span)/ pow(mesh.dy, 2);
//    //exch(af::span, -1, af::span, af::span)+=state.m(af::span, -1, af::span,
//    af::span)/ pow(mesh.dy, 2);
//    //
//    //exch(af::span, af::span, 0 , af::span)+=state.m(af::span, af::span, 0 ,
//    af::span)/ pow(mesh.dz, 2);
//    //exch(af::span, af::span, -1, af::span)+=state.m(af::span, af::span, -1,
//    af::span)/ pow(mesh.dz, 2); if(state.material.afsync) sync();
//    //time_edges += af::timer::stop(timer_edges);
//    computation_time_heff += af::timer::stop(timer_exchsolve);
//    return  (2.* material.A)/(constants::mu0*state.Ms) * exch;
//  }
//  else{
//    timer_exchsolve = af::timer::start();
//    array exch = matmul(matr, flat(state.m));
//    exch=moddims(exch, mesh.n0, mesh.n1, mesh.n2, 3);
//
//    exch.eval();
//    if(state.material.afsync) sync();
//    computation_time_heff += af::timer::stop(timer_exchsolve);
//
//    return  (2.* material.A)/(constants::mu0*state.Ms) * exch;
//  }
//}

// Version yielding real Exchange Field with corrected edges and switch
// conv/sparseMatMul #include "exch.hpp"

//
// void showdims2(const array& a){
//  std::cout<<"Exchange matrix:
//  dims="<<a.dims(0)<<"\t"<<a.dims(1)<<"\t"<<a.dims(2)<<"\t"<<a.dims(3)<<std::endl;
//}
//
////Energy calculation
////Eex=-mu0/2 integral(M . Hex) dx
// double ExchangeField::E(const State& state){
//  return -constants::mu0/2. * state.Ms *
//  afvalue(af::sum(af::sum(af::sum(af::sum(h(state)*state.m, 0), 1), 2), 3)) *
//  mesh.dx * mesh.dy * mesh.dz;
//}
//
////Function returns index
// int ExchangeField::findex(int i0, int i1, int i2, int im, int id){
//  return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
//}
//
////Inner index (index per matrix column)
// int ExchangeField::findex(int i0, int i1, int i2, int im){
//  return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*im));
//}
//
// ExchangeField::ExchangeField (Mesh meshin, Material paramin) :
// material(paramin), mesh(meshin){
//  if(mesh.n0*mesh.n1*mesh.n2>8128){
//    //initialize filters
//    filtr=constant(0.0, 3, 3, 3, f64);
//    filtr(1, 1, 1)= -6 / (pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//
//    filtr(0, 1, 1)= 1 / pow(mesh.dx, 2);
//    filtr(2, 1, 1)= 1 / pow(mesh.dx, 2);
//
//    filtr(1, 0, 1)= 1 / pow(mesh.dy, 2);
//    filtr(1, 2, 1)= 1 / pow(mesh.dy, 2);
//
//    filtr(1, 1, 0)= 1 / pow(mesh.dz, 2);
//    filtr(1, 1, 2)= 1 / pow(mesh.dz, 2);
//  }
//  //Currently better performance for small systems with matmul
//  else{
//    const int dimension=mesh.n0*mesh.n1*mesh.n2*3;
//    double* vmatr = NULL;
//    vmatr = new double[dimension*dimension];
//    for (int id = 0; id < dimension; id++){
//      for (int im = 0; im < 3; im++){
//        for (int i2 = 0; i2 < mesh.n2; i2++){
//          for (int i1 = 0; i1 < mesh.n1; i1++){
//            for (int i0 = 0; i0 < mesh.n0; i0++){
//              const int index=i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
//              vmatr[index]=0.;
//            }
//          }
//        }
//      }
//    }
//
//    for (int id = 0; id < dimension; id++){
//      for (int im = 0; im < 3; im++){
//        for (int i2 = 0; i2 < mesh.n2; i2++){
//          for (int i1 = 0; i1 < mesh.n1; i1++){
//            for (int i0 = 0; i0 < mesh.n0; i0++){
//              const int ind=findex(i0, i1, i2, im);
//              if(ind==id) {
//                vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx,
//                2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
//                //x
//                if(i0==0){
//                  vmatr[findex(i0  , i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                  if (mesh.n0>1) vmatr[findex(i0+1, i1, i2, im,
//                  id)]+= 1./pow(mesh.dx, 2);
//                }
//                if (i0==mesh.n0-1){
//                  vmatr[findex(i0  , i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                  if (mesh.n0>1) vmatr[findex(i0-1, i1, i2, im,
//                  id)]+= 1./pow(mesh.dx, 2);
//                }
//                if(i0>0 && i0< mesh.n0-1){
//                  vmatr[findex(i0-1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                  vmatr[findex(i0+1, i1, i2, im, id)]+= 1./pow(mesh.dx, 2);
//                }
//
//                //y
//                if(i1==0){
//                  vmatr[findex(i0, i1  , i2, im, id)]+= 1./pow(mesh.dy, 2);
//                  if (mesh.n1>1) vmatr[findex(i0, i1+1, i2, im,
//                  id)]+= 1./pow(mesh.dy, 2);
//                }
//                if (i1==mesh.n1-1){
//                  vmatr[findex(i0, i1  , i2, im, id)]+= 1./pow(mesh.dy, 2);
//                  if (mesh.n1>1) vmatr[findex(i0, i1-1, i2, im,
//                  id)]+= 1./pow(mesh.dy, 2);
//                }
//                if(i1>0 && i1< mesh.n1-1){
//                  vmatr[findex(i0, i1-1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//                  vmatr[findex(i0, i1+1, i2, im, id)]+= 1./pow(mesh.dy, 2);
//                }
//
//                //z
//                if (i2==0){
//                  vmatr[findex(i0, i1, i2  , im, id)]+= 1./pow(mesh.dz, 2);
//                  if (mesh.n2>1) vmatr[findex(i0, i1, i2+1, im,
//                  id)]+= 1./pow(mesh.dz, 2);
//                }
//                if (i2==mesh.n2-1){
//                  vmatr[findex(i0, i1, i2  , im, id)]+= 1./pow(mesh.dz, 2);
//                  if (mesh.n2>1) vmatr[findex(i0, i1, i2-1, im,
//                  id)]+= 1./pow(mesh.dz, 2);
//                }
//                if(i2>0 && i2< mesh.n2-1){
//                  vmatr[findex(i0, i1, i2-1, im, id)]+= 1./pow(mesh.dz, 2);
//                  vmatr[findex(i0, i1, i2+1, im, id)]+= 1./pow(mesh.dz, 2);
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//    array fullmatr(dimension, dimension, vmatr);
//    delete [] vmatr;
//    vmatr = NULL;
//    //showdims2(fullmatr);
//    matr=sparse(fullmatr);
//
//    std::cout << "Sparsity of matr = "
//              << (float)sparseGetNNZ(matr) / (float)matr.elements()
//              << std::endl;
//    }
//}
//
// array ExchangeField::h(const State& state){
//  timer_exchsolve = af::timer::start();
//
//  if(mesh.n0*mesh.n1*mesh.n2>8128){
//    timer_conv = af::timer::start();
//    //convolution
//    array exch = convolve(state.m, filtr, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
//
//    if(state.material.afsync) sync();
//    time_conv += af::timer::stop(timer_conv);
//
//    //Accounting for boundary conditions by adding initial m values on the
//    boundaries by adding all 6 boundary surfaces timer_edges =
//    af::timer::start(); exch(0, af::span, af::span, af::span)+=state.m(0 ,
//    af::span, af::span, af::span)/ pow(mesh.dx, 2); exch(-1, af::span,
//    af::span, af::span)+=state.m(-1, af::span, af::span, af::span)/
//    pow(mesh.dx, 2);
//
//
//    exch(af::span, 0 , af::span, af::span)+=state.m(af::span, 0 , af::span,
//    af::span)/ pow(mesh.dy, 2); exch(af::span, -1, af::span,
//    af::span)+=state.m(af::span, -1, af::span, af::span)/ pow(mesh.dy, 2);
//
//    exch(af::span, af::span, 0 , af::span)+=state.m(af::span, af::span, 0 ,
//    af::span)/ pow(mesh.dz, 2); exch(af::span, af::span, -1,
//    af::span)+=state.m(af::span, af::span, -1, af::span)/ pow(mesh.dz, 2);
//    if(state.material.afsync) sync();
//    time_edges += af::timer::stop(timer_edges);
//    computation_time_heff += af::timer::stop(timer_exchsolve);
//    return  (2.* material.A)/(constants::mu0*state.Ms) * exch;
//  }
//  else{
//    timer_exchsolve = af::timer::start();
//    array exch = matmul(matr, flat(state.m));
//    exch=moddims(exch, mesh.n0, mesh.n1, mesh.n2, 3);
//
//    exch.eval();
//    if(state.material.afsync) sync();
//    computation_time_heff += af::timer::stop(timer_exchsolve);
//
//    return  (2.* material.A)/(constants::mu0*state.Ms) * exch;
//  }
//}
} // namespace magnumafcpp
