#include "micro_exch_sparse.hpp"
#include "../misc.hpp"
#include "../func.hpp"

namespace magnumaf{


SparseExchangeField::SparseExchangeField (float A_exchange, Mesh mesh, bool verbose) : matr(calc_CSR_matrix(A_exchange, mesh, verbose))
{
}


SparseExchangeField::SparseExchangeField (const af::array& A_exchange_field, Mesh mesh, bool verbose) : matr(calc_CSR_matrix(A_exchange_field, mesh, verbose))
{
}


// For wrapping only: constructor version taking A_exchange_field
SparseExchangeField::SparseExchangeField (long int A_exchange_field_ptr, Mesh mesh, bool verbose) : matr(calc_CSR_matrix(*(new af::array( *((void **) A_exchange_field_ptr))), mesh, verbose))
{
}


af::array SparseExchangeField::h(const State& state){
    af::timer aftimer = af::timer::start();
    af::array exch = af::matmul(matr, af::flat(state.m));
    exch = af::moddims(exch, state.mesh.n0, state.mesh.n1, state.mesh.n2, 3);
    if(state.afsync) af::sync();
    af_time += af::timer::stop(aftimer);
    if (state.Ms_field.isempty()){
        return  exch/state.Ms;
    }
    else {
        af::array heff = exch/state.Ms_field;
        replace(heff, state.Ms_field!=0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}


// Energy calculation: E_ex = -mu0/2 * integral(M * Hex) dx
float SparseExchangeField::E(const State& state){
    if( state.Ms_field.isempty() ){
        return -constants::mu0/2. * state.Ms * afvalue(af::sum(af::sum(af::sum(af::sum( h(state) * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
    }
    else{
        return -constants::mu0/2. * afvalue(af::sum(af::sum(af::sum(af::sum(state.Ms_field * h(state) * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
    }
}


float SparseExchangeField::E(const State& state, const af::array& h){
    if( state.Ms_field.isempty() ){
        return -constants::mu0/2. * state.Ms * afvalue(af::sum(af::sum(af::sum(af::sum( h * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
    }
    else{
        return -constants::mu0/2. * afvalue(af::sum(af::sum(af::sum(af::sum(state.Ms_field * h * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
    }
}


// Get inner index (index per matrix column)
int SparseExchangeField::findex(int i0, int i1, int i2, int im, Mesh mesh){
    return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*im));
}


af::array SparseExchangeField::calc_CSR_matrix(const float A_exchange, const Mesh& mesh, const bool verbose){
    af::timer t;
    if(verbose) af::timer::start();
    const int dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;

    std::vector<float> CSR_values;// matrix values,  of length "number of elements"
    std::vector<int> CSR_IA (dimension + 1);// recursive row indices of length (n_rows + 1): IA[0] = 0; IA[i] = IA[i-1] + (number of nonzero elements on the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;// comumn index of each element, hence of length "number of elements"
    for (int id = 0; id < dimension; id++){// loop over rows (or cols?^^)
      int csr_ia = 0; // counter for SCR_IA
      for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.n2; i2++){
          for (int i1 = 0; i1 < mesh.n1; i1++){
            for (int i0 = 0; i0 < mesh.n0; i0++){
              const int ind=findex(i0, i1, i2, im, mesh);
              if(ind==id) {
                //std::cout << ind << ", " << id << ", " << im << ", " << i2 << ", " << i1 << ", " << i0 << std::endl;
                //Note: skippable due to cross product property://vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
                //x
                if(i0 == 0 && mesh.n0 > 1 ){
                    CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow( mesh.dx, 2) );
                    CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                    csr_ia++;
                }
                if (i0 == mesh.n0 - 1 && mesh.n0 > 1){
                    CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dx, 2) );
                    CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                    csr_ia++;
                }
                if(i0>0 && i0< mesh.n0 - 1 ){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dx, 2) );
                  CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                  csr_ia++;

                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dx, 2) );
                  CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                  csr_ia++;
                }

                //y
                if(i1 == 0 && mesh.n1 > 1 ){
                    CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                    CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                    csr_ia++;
                }
                if (i1 == mesh.n1 - 1 && mesh.n1 > 1){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                  CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                  csr_ia++;
                }
                if(i1>0 && i1< mesh.n1-1){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                  CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                  csr_ia++;
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                  CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                  csr_ia++;
                }

                //z
                if (i2 == 0 && mesh.n2 > 1 ){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dz, 2) );
                  CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                  csr_ia++;
                }
                if (i2 == mesh.n2 - 1 && mesh.n2 > 1){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dz, 2) );
                  CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                  csr_ia++;
                }
                if( i2 > 0 && i2 < mesh.n2 - 1){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dz, 2) );
                  CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                  csr_ia++;
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dz, 2) );
                  CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                  csr_ia++;
                }
              }
            }
          }
        }
      }
      CSR_IA[id + 1] = CSR_IA[id] + csr_ia;
    }

    af::array result = af::sparse((dim_t) dimension, (dim_t) dimension, (dim_t) CSR_values.size(), (void*) CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f32);
    if(verbose) printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of CSR_matrix = %f\n", Info(), t.stop(), (float)af::sparseGetNNZ(matr) / (float)matr.elements());
    return result;
}

// Assembly of sparse matrix for spacially varying exchange energy A_exchange_field
af::array SparseExchangeField::calc_CSR_matrix(const af::array& A_exchange_field, const Mesh& mesh, const bool verbose){
    printf("%s SparseExchangeField::calc_CSR_matrix unit testing not finished!\n", Warning());
    fflush(stdout);
    af::timer t;
    if(verbose) af::timer::start();
    const int dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;

    std::vector<float> CSR_values;// matrix values,  of length "number of elements"
    std::vector<int> CSR_IA (dimension + 1);// recursive row indices of length (n_rows + 1): IA[0] = 0; IA[i] = IA[i-1] + (number of nonzero elements on the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;// comumn index of each element, hence of length "number of elements"
    float* a_host = NULL;
    a_host = A_exchange_field.host<float>();
    for (int id = 0; id < dimension; id++){// loop over rows (or cols?^^)
      int csr_ia = 0; // counter for SCR_IA
      for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.n2; i2++){
          for (int i1 = 0; i1 < mesh.n1; i1++){
            for (int i0 = 0; i0 < mesh.n0; i0++){
              const int ind=findex(i0, i1, i2, im, mesh);
              if(ind==id) {
                //std::cout << ind << ", " << id << ", " << im << ", " << i2 << ", " << i1 << ", " << i0 << std::endl;
                //Note: skippable due to cross product property://vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
                // x
                // Note: poor indexing performace. TODO improve performance: directly accessing values with afvalue increades sp4 assembly from ~0.4 s to ~1.4 s! maybe access full host array once?
                // is host data then in correct order for adapted findex for scalar field, i.e. i0 + mesh.n0 * (i1 + mesh.n1 * i2)?
                // TODO consider changing A_exchange_field(i0+1, i1, i2) to 'local' A_exchange_field(i0, i1, i2) for x, y, z
                if(i0 == 0 && mesh.n0 > 1){
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_p = a_host[util::stride(i0+1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                        csr_ia++;
                    }
                }
                else if (i0 == mesh.n0 - 1 && mesh.n0 > 1){
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_m = a_host[util::stride(i0-1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if(i0>0 && i0< mesh.n0 - 1 ){
                    //i_x +u1
                    {
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_m = a_host[util::stride(i0-1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i_m != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                        csr_ia++;
                    }
                    }

                    {
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_p = a_host[util::stride(i0+1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i_p != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                        csr_ia++;
                    }
                    }
                }

                // y
                if(i1 == 0 && mesh.n1 > 1 ){
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_p = a_host[util::stride(i0, i1+1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if (i1 == mesh.n1 - 1 && mesh.n1 > 1){
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_m = a_host[util::stride(i0, i1-1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if(i1>0 && i1< mesh.n1-1){
                    {
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_m = a_host[util::stride(i0, i1-1, i2, mesh.n0, mesh.n1)];
                    if (A_i_m != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                        csr_ia++;
                    }
                    }
                    {
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_p = a_host[util::stride(i0, i1+1, i2, mesh.n0, mesh.n1)];
                    if (A_i_p != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                        csr_ia++;
                    }
                    }
                }

                // z
                if (i2 == 0 && mesh.n2 > 1 ){
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_p = a_host[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if (i2 == mesh.n2 - 1 && mesh.n2 > 1){
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_m = a_host[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if( i2 > 0 && i2 < mesh.n2 - 1){
                    {
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_m = a_host[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];
                    if (A_i_m != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                        csr_ia++;
                    }
                    }
                    {
                    float A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    float A_i_p = a_host[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                    if (A_i_p != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                        csr_ia++;
                    }
                    }
                }
              }
            }
          }
        }
      }
      CSR_IA[id + 1] = CSR_IA[id] + csr_ia;
    }
    af::freeHost(a_host);
    af::array result = af::sparse((dim_t) dimension, (dim_t) dimension, (dim_t) CSR_values.size(), (void*) CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f32);
    if(verbose) {
        printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of CSR_matrix = %f\n", Info(), t.stop(), (float)af::sparseGetNNZ(result) / (float)result.elements());
        fflush(stdout);
    }
    return result;
}
}// namespace magnumaf
