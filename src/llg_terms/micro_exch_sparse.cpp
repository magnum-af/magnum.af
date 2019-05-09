#include "micro_exch_sparse.hpp"

//Energy calculation
//Eex=-mu0/2 integral(M . Hex) dx
double SparseExchangeField::E(const State& state){
    return -constants::mu0/2. * state.material.ms * afvalue(af::sum(af::sum(af::sum(af::sum(h(state)*state.m,0),1),2),3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz; 
}

double SparseExchangeField::E(const State& state, const af::array& h){
    return -constants::mu0/2. * state.material.ms * afvalue(sum(sum(sum(sum(h * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}

//Function returns index 
int SparseExchangeField::findex(int i0, int i1, int i2, int im, int id, Mesh mesh){
    return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
}

//Inner index (index per matrix column)
int SparseExchangeField::findex(int i0, int i1, int i2, int im, Mesh mesh){
    return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*im));
}

SparseExchangeField::SparseExchangeField (Mesh mesh, Material material){
    const int dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;
    double* vmatr = NULL;
    vmatr = new double[dimension*dimension];
    for (int id = 0; id < dimension; id++){
      for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.n2; i2++){
          for (int i1 = 0; i1 < mesh.n1; i1++){
            for (int i0 = 0; i0 < mesh.n0; i0++){
              const int index=i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*(im+3*id)));
              vmatr[index]=0.;
            }
          }
        }
      }
    }
  
    for (int id = 0; id < dimension; id++){
      for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.n2; i2++){
          for (int i1 = 0; i1 < mesh.n1; i1++){
            for (int i0 = 0; i0 < mesh.n0; i0++){
              const int ind=findex(i0,i1,i2,im, mesh);
              if(ind==id) {
                //Note: skippable due to cross product property://vmatr[findex(i0,i1,i2,im,id)]+=-6./(pow(mesh.dx,2)+pow(mesh.dy,2)+pow(mesh.dz,2));
                //x
                if(i0==0){
                  //skippable//vmatr[findex(i0  ,i1,i2,im,id)]+= 1./pow(mesh.dx,2);
                  if (mesh.n0>1) vmatr[findex(i0+1,i1,i2,im,id, mesh)]+= 1./pow(mesh.dx,2);
                }
                if (i0==mesh.n0-1){
                  //skippable//vmatr[findex(i0  ,i1,i2,im,id)]+= 1./pow(mesh.dx,2);
                  if (mesh.n0>1) vmatr[findex(i0-1,i1,i2,im,id, mesh)]+= 1./pow(mesh.dx,2);
                }
                if(i0>0 && i0< mesh.n0-1){
                  vmatr[findex(i0-1,i1,i2,im,id, mesh)]+= 1./pow(mesh.dx,2);
                  vmatr[findex(i0+1,i1,i2,im,id, mesh)]+= 1./pow(mesh.dx,2);
                }
  
                //y
                if(i1==0){
                  //skippable//vmatr[findex(i0,i1  ,i2,im,id)]+= 1./pow(mesh.dy,2);
                  if (mesh.n1>1) vmatr[findex(i0,i1+1,i2,im,id, mesh)]+= 1./pow(mesh.dy,2);
                }
                if (i1==mesh.n1-1){
                  //skippable//vmatr[findex(i0,i1  ,i2,im,id)]+= 1./pow(mesh.dy,2);
                  if (mesh.n1>1) vmatr[findex(i0,i1-1,i2,im,id, mesh)]+= 1./pow(mesh.dy,2);
                }                     
                if(i1>0 && i1< mesh.n1-1){
                  vmatr[findex(i0,i1-1,i2,im,id, mesh)]+= 1./pow(mesh.dy,2);
                  vmatr[findex(i0,i1+1,i2,im,id, mesh)]+= 1./pow(mesh.dy,2);
                }
  
                //z
                if (i2==0){
                  //skippable//vmatr[findex(i0,i1,i2  ,im,id)]+= 1./pow(mesh.dz,2);
                  if (mesh.n2>1) vmatr[findex(i0,i1,i2+1,im,id, mesh)]+= 1./pow(mesh.dz,2);
                }
                if (i2==mesh.n2-1){
                  //skippable//vmatr[findex(i0,i1,i2  ,im,id)]+= 1./pow(mesh.dz,2);
                  if (mesh.n2>1) vmatr[findex(i0,i1,i2-1,im,id, mesh)]+= 1./pow(mesh.dz,2);
                }
                if(i2>0 && i2< mesh.n2-1){
                  vmatr[findex(i0,i1,i2-1,im,id, mesh)]+= 1./pow(mesh.dz,2);
                  vmatr[findex(i0,i1,i2+1,im,id, mesh)]+= 1./pow(mesh.dz,2);
                }
              }
            }
          }
        }
      }
    }
    af::array fullmatr(dimension, dimension, vmatr);
    delete [] vmatr;
    vmatr = NULL;
    matr = af::sparse(fullmatr);
  
    std::cout << "Sparsity of matr = "
              << (float)af::sparseGetNNZ(matr) / (float)matr.elements()
              << std::endl;
}

af::array SparseExchangeField::h(const State& state){
    af::timer aftimer = af::timer::start();
    af::array exch = af::matmul(matr, af::flat(state.m));
    exch = moddims(exch, state.mesh.n0, state.mesh.n1, state.mesh.n2, 3);
    exch.eval();//TODO check if can be skipped
    if(state.material.afsync) sync();
    af_time += af::timer::stop(aftimer);
    return  (2.* state.material.A)/(constants::mu0 * state.material.ms) * exch;
}
