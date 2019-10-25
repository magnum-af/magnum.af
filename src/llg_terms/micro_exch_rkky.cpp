#include "micro_exch_rkky.hpp"
#include "../misc.hpp"
#include "../func.hpp"

namespace magnumafcpp{


RKKYExchangeField::RKKYExchangeField (RKKY_values rkky_values, Exchange_values exchange_values, Mesh mesh, const af::array& rkky_indices, bool verbose) : matr(calc_CSR_matrix(rkky_values.get(), exchange_values.get(), mesh, rkky_indices, verbose))
{
}


RKKYExchangeField::RKKYExchangeField (int todelswitch, RKKY_values rkky_values, Exchange_values exchange_values, Mesh mesh, const af::array& rkky_indices) : matr(calc_COO_matrix(rkky_values.get(), exchange_values.get(), mesh, rkky_indices, true))
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
af::array RKKYExchangeField::calc_CSR_matrix(const af::array& RKKY_field, const af::array& A_exchange_field, const Mesh& mesh, const af::array& rkky_indices, const bool verbose){
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
    unsigned int* rkky_indices_raw = NULL;
    //af::print("", rkky_indices);
    //rkky_indices_raw = rkky_indices.host<unsigned int>();
    if ( ! rkky_indices.isempty() ) rkky_indices_raw = rkky_indices.host<unsigned int>();
    for (int id = 0; id < dimension; id++){// loop over rows (or cols?^^)
      int csr_ia = 0; // counter for SCR_IA
      //openmp for docker before//
      //TODO//#pragma omp parallel for
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

                    //const unsigned int RKKY_index_i = rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    //const unsigned int RKKY_index_i_p = rkky_indices_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i_p = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                    //std::cout << "rkkyindex = " << RKKY_index_i << "and" << RKKY_index_i_p << std::endl;
                    if ( (RKKY_index_i == RKKY_index_i_p) && (RKKY_i != 0) && ( RKKY_i_p != 0) ){
                        //assuming rkky jump condition equal to exch jump
                        CSR_values.push_back( (2.* RKKY_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* RKKY_i_p/(RKKY_i_p + RKKY_i));
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

                    //const unsigned int RKKY_index_i = rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    //const unsigned int RKKY_index_i_m = rkky_indices_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i_m = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];

                    if ( (RKKY_index_i == RKKY_index_i_m) && (RKKY_i != 0) && ( RKKY_i_m != 0)){
                        //assuming rkky jump condition equal to exch jump
                        CSR_values.push_back( (2.* RKKY_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* RKKY_i_m/(RKKY_i_m + RKKY_i));
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
    //for (auto const& value: CSR_IA){
    //    std::cout << "CSR_IA=" << value << std::endl;
    //}
    //for (auto const& value: CSR_JA){
    //    std::cout << "CSR_JA=" << value << std::endl;
    //}
    af::freeHost(a_raw);
    af::freeHost(rkky_raw);
    af::freeHost(rkky_indices_raw);
    af::array result = af::sparse((dim_t) dimension, (dim_t) dimension, (dim_t) CSR_values.size(), (void*) CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if(verbose) {
        printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of CSR_matrix = %f\n", Info(), t.stop(), (float)af::sparseGetNNZ(result) / (float)result.elements());
        fflush(stdout);
    }
    return result;
}


// Assembly of sparse matrix for spacially varying exchange energy A_exchange_field
af::array RKKYExchangeField::calc_COO_matrix(const af::array& RKKY_field, const af::array& A_exchange_field, const Mesh& mesh, const af::array& rkky_indices, const bool verbose){
    printf("%s RKKYExchangeField::calc_COO_matrix unit testing not finished!\n", Warning());
    fflush(stdout);
    af::timer t;
    if(verbose) af::timer::start();
    const int dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;
    std::vector<double> COO_values;// matrix values,  of length "number of elements"
    std::vector<int> COO_COL;
    std::vector<int> COO_ROW;
    double* a_raw = NULL;
    a_raw = A_exchange_field.host<double>();
    double* rkky_raw = NULL;
    rkky_raw = RKKY_field.host<double>();
    unsigned int* rkky_indices_raw = NULL;
    //af::print("", rkky_indices);
    //rkky_indices_raw = rkky_indices.host<unsigned int>();
    if ( ! rkky_indices.isempty() ) rkky_indices_raw = rkky_indices.host<unsigned int>();
    //openmp for docker before//
    //TODO aborts program//#pragma omp parallel for
    for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.n2; i2++){
            for (int i1 = 0; i1 < mesh.n1; i1++){
                for (int i0 = 0; i0 < mesh.n0; i0++){
                    const int ind=findex(i0, i1, i2, im, mesh);
                    //TODO check if COL and ROW are correct ordered
                // x
                if( (i0 == 0 && mesh.n0 > 1) || (i0>0 && i0< mesh.n0 - 1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_p = a_raw[util::stride(i0+1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        COO_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p/(A_i_p + A_i));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back( findex( i0+1, i1, i2, im, mesh) );
                    }
                }
                if ( (i0 == mesh.n0 - 1 && mesh.n0 > 1) || (i0>0 && i0< mesh.n0 - 1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_m = a_raw[util::stride(i0-1, i1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        COO_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back( findex( i0-1, i1, i2, im, mesh ) );
                    }
                }

                // y
                if( (i1 == 0 && mesh.n1 > 1) || (i1>0 && i1< mesh.n1-1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_p = a_raw[util::stride(i0, i1+1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        COO_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back( findex( i0, i1+1, i2, im, mesh ) );
                    }
                }
                if ( (i1 == mesh.n1 - 1 && mesh.n1 > 1) || (i1>0 && i1< mesh.n1-1)){
                    double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double A_i_m = a_raw[util::stride(i0, i1-1, i2, mesh.n0, mesh.n1)];
                    if (A_i != 0){
                        COO_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back( findex( i0, i1-1, i2, im, mesh ) );
                    }
                }

                // z
                // Preferring RKKY over exch vals
                if ( (i2 == 0 && mesh.n2 > 1) || ( i2 > 0 && i2 < mesh.n2 - 1)){
                    double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double RKKY_i_p = rkky_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];

                    //const unsigned int RKKY_index_i = rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    //const unsigned int RKKY_index_i_p = rkky_indices_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i_p = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                    //std::cout << "rkkyindex = " << RKKY_index_i << "and" << RKKY_index_i_p << std::endl;
                    if ( (RKKY_index_i == RKKY_index_i_p) && (RKKY_i != 0) && ( RKKY_i_p != 0) ){
                        //assuming rkky jump condition equal to exch jump
                        COO_values.push_back( (2.* RKKY_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* RKKY_i_p/(RKKY_i_p + RKKY_i));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back( findex(i0, i1, i2+1, im, mesh) );
                    }
                    else{
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_raw[util::stride(i0, i1, i2+1, mesh.n0, mesh.n1)];
                        if (A_i != 0){
                            COO_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_p/(A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back( findex( i0, i1, i2+1, im, mesh ) );
                        }
                    }
                }
                if ( (i2 == mesh.n2 - 1 && mesh.n2 > 1) || ( i2 > 0 && i2 < mesh.n2 - 1)){
                    double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    double RKKY_i_m = rkky_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];

                    //const unsigned int RKKY_index_i = rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    //const unsigned int RKKY_index_i_m = rkky_indices_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                    const unsigned int RKKY_index_i_m = rkky_indices.isempty()? 0 : rkky_indices_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];

                    if ( (RKKY_index_i == RKKY_index_i_m) && (RKKY_i != 0) && ( RKKY_i_m != 0)){
                        //assuming rkky jump condition equal to exch jump
                        COO_values.push_back( (2.* RKKY_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* RKKY_i_m/(RKKY_i_m + RKKY_i));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back( findex( i0, i1, i2-1, im, mesh ) );
                    }
                    else{
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_raw[util::stride(i0, i1, i2-1, mesh.n0, mesh.n1)];
                        if (A_i != 0){
                            COO_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dz, 2)) * 2.* A_i_m/(A_i_m + A_i));
                                COO_ROW.push_back(ind);
                                COO_COL.push_back( findex( i0, i1, i2-1, im, mesh ) );
                        }
                    }
                }
                }
            }
        }
    }

    ////TODO remove
    //std::cout << "test" << std::endl;
    //COO_values.push_back(1.);
    //COO_COL.push_back(0);
    //COO_ROW.push_back(0);
    //af::print("", af::array(COO_values.size(), COO_values.data()));
    ////causes error //af::print("", af::array(COO_values.data()));
    ////end TODO remove

    af::freeHost(a_raw);
    af::freeHost(rkky_raw);
    af::freeHost(rkky_indices_raw);
    af::array result = af::sparse((dim_t) dimension, (dim_t) dimension, af::array(COO_values.size(), COO_values.data()), af::array(COO_ROW.size(), COO_ROW.data()), af::array(COO_COL.size(), COO_COL.data()), AF_STORAGE_COO);
    double time = t.stop();

    af::timer timer_convert = af::timer::start();
    result = af::sparseConvertTo(result, AF_STORAGE_CSR);
    double time_convert = timer_convert.stop();

    if(verbose) {
        printf("%s Initialized sparse RKKY-exchange matrix in %f [s]. Converted COO to CSR in %f [s]. Sparsity = %f\n", Info(), time, time_convert, (float)af::sparseGetNNZ(result) / (float)result.elements());
        fflush(stdout);
    }
    //moved up for timing//return af::sparseConvertTo(result, AF_STORAGE_CSR);
    return result;
}


}// namespace magnumafcpp
