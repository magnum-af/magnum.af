#include "micro_exch_rkky.hpp"
#include "../misc.hpp"
#include "../func.hpp"

namespace magnumaf{


RKKYExchangeField::RKKYExchangeField (RKKY_values rkky_values, Exchange_values exchange_values, Mesh mesh, bool verbose) : matr(calc_CSR_matrix(rkky_values.get(), exchange_values.get(), mesh, verbose))
{
}


af::array RKKYExchangeField::h(const State& state){
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


// Get inner index (index per matrix column)
int RKKYExchangeField::findex(int i0, int i1, int i2, int im, Mesh mesh){
    return i0+mesh.n0*(i1+mesh.n1*(i2+mesh.n2*im));
}


// Assembly of sparse matrix for spacially varying exchange energy A_exchange_field
af::array RKKYExchangeField::calc_CSR_matrix(const af::array& RKKY_field, const af::array& A_exchange_field, const Mesh& mesh, const bool verbose){
    printf("%s RKKYExchangeField::calc_CSR_matrix unit testing not finished!\n", Warning());
    fflush(stdout);
    af::timer t;
    if(verbose) af::timer::start();
    const int dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;

    std::vector<double> CSR_values;// matrix values,  of length "number of elements"
    std::vector<int> CSR_IA (dimension + 1);// recursive row indices of length (n_rows + 1): IA[0] = 0; IA[i] = IA[i-1] + (number of nonzero elements on the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;// comumn index of each element, hence of length "number of elements"
    double* a_raw = NULL;
    a_raw = A_exchange_field.host<double>();
    double* rkky_raw = NULL;
    rkky_raw = RKKY_field.host<double>();
    for (int id = 0; id < dimension; id++){// loop over rows (or cols?^^)
      int csr_ia = 0; // counter for SCR_IA
      for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.n2; i2++){
          for (int i1 = 0; i1 < mesh.n1; i1++){
            for (int i0 = 0; i0 < mesh.n0; i0++){
              const int ind=findex(i0, i1, i2, im, mesh);
              if(ind==id) {
                //Note: skippable due to cross product property://vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
                // x
                if( (i0 == 0 && mesh.n0 > 1) || (i0>0 && i0< mesh.n0 - 1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_p = a_raw[util::stride(i0+1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                        csr_ia++;
                    }
                }
                if ( (i0 == mesh.n0 - 1 && mesh.n0 > 1) || (i0>0 && i0< mesh.n0 - 1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_m = a_raw[util::stride(i0-1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }

                // y
                if( (i1 == 0 && mesh.n1 > 1) || (i1>0 && i1< mesh.n1-1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_p = a_raw[util::stride(i0, i1+1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }
                if ( (i1 == mesh.n1 - 1 && mesh.n1 > 1) || (i1>0 && i1< mesh.n1-1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_m = a_raw[util::stride(i0, i1-1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }

                // z
                // Preferring RKKY over exch vals
                if ( (i2 == 0 && mesh.n2 > 1) || ( i2 > 0 && i2 < mesh.n2 - 1)){
                    double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double RKKY_i_p = rkky_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                    if ( (RKKY_i != 0) && ( RKKY_i_p != 0)){
                        CSR_values.push_back( (RKKY_i + RKKY_i_p)/2. );//NOTE: maybe other norm?
                        CSR_JA.push_back( findex(i0, i1, i2+1, im, mesh) );
                        csr_ia++;
                    }
                    else{
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                        if (A_i != 0){
                            CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_p/(A_i_p + A_i));
                            CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                            csr_ia++;
                        }
                    }
                }
                if ( (i2 == mesh.n2 - 1 && mesh.n2 > 1) || ( i2 > 0 && i2 < mesh.n2 - 1)){
                    double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double RKKY_i_m = rkky_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];

                    if ( (RKKY_i != 0) && ( RKKY_i_m != 0)){
                        CSR_values.push_back( (RKKY_i + RKKY_i_m)/2. );
                        CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                        csr_ia++;
                    }
                    else{
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];
                        if (A_i != 0){
                            CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_m/(A_i_m + A_i));
                            CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
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
    af::freeHost(a_raw);
    af::freeHost(rkky_raw);
    af::array result = af::sparse((dim_t) dimension, (dim_t) dimension, (dim_t) CSR_values.size(), (void*) CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if(verbose) {
        printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of CSR_matrix = %f\n", Info(), t.stop(), (float)af::sparseGetNNZ(result) / (float)result.elements());
        fflush(stdout);
    }
    return result;
}
}// namespace magnumaf
