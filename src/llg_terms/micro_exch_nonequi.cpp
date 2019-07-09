#include "micro_exch_nonequi.hpp"
#include "../misc.hpp"
#include "../func.hpp"

namespace magnumaf{


NonequiExchangeField::NonequiExchangeField (double A_exchange, NonequispacedMesh mesh, bool verbose) : matr(calc_CSR_matrix(A_exchange, mesh, verbose))
{
}


NonequiExchangeField::NonequiExchangeField (const af::array& A_exchange_field, NonequispacedMesh mesh, bool verbose) : matr(calc_CSR_matrix(A_exchange_field, mesh, verbose))
{
}


// For wrapping only: constructor version taking A_exchange_field
NonequiExchangeField::NonequiExchangeField (long int A_exchange_field_ptr, NonequispacedMesh mesh, bool verbose) : matr(calc_CSR_matrix(*(new af::array( *((void **) A_exchange_field_ptr))), mesh, verbose))
{
}


af::array NonequiExchangeField::h(const State& state){
    af::timer aftimer = af::timer::start();
    af::array exch = af::matmul(matr, af::flat(state.m));
    exch = af::moddims(exch, state.nonequimesh.nx, state.nonequimesh.ny, state.nonequimesh.nz, 3);
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
double NonequiExchangeField::E(const State& state){
    return -constants::mu0/2. * state.Ms * afvalue(af::sum(af::sum(af::sum(af::sum(h(state)*state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}


double NonequiExchangeField::E(const State& state, const af::array& h){
    return -constants::mu0/2. * state.Ms * afvalue(sum(sum(sum(sum(h * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
}
// Get inner index (index per matrix column)
int NonequiExchangeField::findex(int i0, int i1, int i2, int im, NonequispacedMesh mesh){
    return i0+mesh.nx*(i1+mesh.ny*(i2+mesh.nz*im));
}


af::array NonequiExchangeField::calc_CSR_matrix(const double A_exchange, const NonequispacedMesh& mesh, const bool verbose){
    printf("%s NonequiExchangeField::calc_CSR_matrix unit testing not finished!\n", Warning()); fflush(stdout);
    af::timer t;
    if(verbose) af::timer::start();

    std::vector<double> h;// spacings between discretization points h = (dz[n] + dz[n+1])/2
    for(unsigned int i = 0; i < mesh.z_spacing.size() - 1; i++){
        h.push_back((mesh.z_spacing.at(i) + mesh.z_spacing.at(i+1))/2.);
    }

    const int dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> CSR_values;// matrix values,  of length "number of elements"
    std::vector<int> CSR_IA (dimension + 1);// recursive row indices of length (n_rows + 1): IA[0] = 0; IA[i] = IA[i-1] + (number of nonzero elements on the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;// comumn index of each element, hence of length "number of elements"
    for (int id = 0; id < dimension; id++){// loop over rows (or cols?^^)
      int csr_ia = 0; // counter for SCR_IA
      for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.nz; i2++){
          for (int i1 = 0; i1 < mesh.ny; i1++){
            for (int i0 = 0; i0 < mesh.nx; i0++){
              const int ind=findex(i0, i1, i2, im, mesh);
              if(ind==id) {
                //std::cout << ind << ", " << id << ", " << im << ", " << i2 << ", " << i1 << ", " << i0 << std::endl;
                //Note: skippable due to cross product property://vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
                //x
                if(i0 == 0 && mesh.nx > 1 ){
                    CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow( mesh.dx, 2) );
                    CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                    csr_ia++;
                }
                if (i0 == mesh.nx - 1 && mesh.nx > 1){
                    CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dx, 2) );
                    CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                    csr_ia++;
                }
                if(i0>0 && i0< mesh.nx - 1 ){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dx, 2) );
                  CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                  csr_ia++;

                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dx, 2) );
                  CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                  csr_ia++;
                }

                //y
                if(i1 == 0 && mesh.ny > 1 ){
                    CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                    CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                    csr_ia++;
                }
                if (i1 == mesh.ny - 1 && mesh.ny > 1){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                  CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                  csr_ia++;
                }
                if(i1>0 && i1< mesh.ny-1){
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                  CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                  csr_ia++;
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(mesh.dy, 2) );
                  CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                  csr_ia++;
                }

                // z from ref [1]
                if (i2 == 0 && mesh.nz > 1 ){//TODO check
                  // Note: skipping f-1 term as it drops out in llg: neumann bc is assumed, which would consider fictive m[-1] with value m[0]
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(h.at(i2), 2) );
                  CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                  csr_ia++;
                }
                else if (i2 == mesh.nz - 1 && mesh.nz > 1){//TODO check
                  // Note: skipping f+1 term as it drops out in llg: neumann bc is assumed, which would consider fictive m[n] with value m[n-1]
                  CSR_values.push_back( (2.* A_exchange)/(constants::mu0) * 1./pow(h.at(i2-1), 2) );
                  CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                  csr_ia++;
                }
                else if( i2 > 0 && i2 < mesh.nz - 1){
                  double h_divisor = h.at(i2) * h.at(i2-1) * (1. + h.at(i2)/h.at(i2-1));
                  // f_{i-1} term
                  CSR_values.push_back( (2.* A_exchange / constants::mu0) * (2. * h.at(i2)/h.at(i2-1))/h_divisor);
                  CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                  csr_ia++;
                  // f_{i+1} term
                  CSR_values.push_back( (2.* A_exchange / constants::mu0) * 2./h_divisor );
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

    af::array result = af::sparse((dim_t) dimension, (dim_t) dimension, (dim_t) CSR_values.size(), (void*) CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if(verbose) printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of CSR_matrix = %f\n", Info(), t.stop(), (float)af::sparseGetNNZ(matr) / (float)matr.elements());
    return result;
}

// Assembly of sparse matrix for spacially varying exchange energy A_exchange_field
af::array NonequiExchangeField::calc_CSR_matrix(const af::array& A_exchange_field, const NonequispacedMesh& mesh, const bool verbose){
    printf("%s NonequiExchangeField::calc_CSR_matrix unit testing not finished!\n", Warning()); fflush(stdout);
    af::timer t;
    if(verbose) af::timer::start();

    std::vector<double> h;// spacings between discretization points h = (dz[n] + dz[n+1])/2
    for(unsigned int i = 0; i < mesh.z_spacing.size() - 1; i++){
        h.push_back((mesh.z_spacing.at(i) + mesh.z_spacing.at(i+1))/2.);
    }

    const int dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> CSR_values;// matrix values,  of length "number of elements"
    std::vector<int> CSR_IA (dimension + 1);// recursive row indices of length (n_rows + 1): IA[0] = 0; IA[i] = IA[i-1] + (number of nonzero elements on the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;// comumn index of each element, hence of length "number of elements"
    double* a_host = NULL;
    a_host = A_exchange_field.host<double>();
    for (int id = 0; id < dimension; id++){// loop over rows (or cols?^^)
      int csr_ia = 0; // counter for SCR_IA
      for (int im = 0; im < 3; im++){
        for (int i2 = 0; i2 < mesh.nz; i2++){
          for (int i1 = 0; i1 < mesh.ny; i1++){
            for (int i0 = 0; i0 < mesh.nx; i0++){
              const int ind=findex(i0, i1, i2, im, mesh);
              if(ind==id) {
                //std::cout << ind << ", " << id << ", " << im << ", " << i2 << ", " << i1 << ", " << i0 << std::endl;
                //Note: skippable due to cross product property://vmatr[findex(i0, i1, i2, im, id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz, 2));
                // x
                // Note: poor indexing performace. TODO improve performance: directly accessing values with afvalue increades sp4 assembly from ~0.4 s to ~1.4 s! maybe access full host array once?
                // is host data then in correct order for adapted findex for scalar field, i.e. i0 + mesh.nx * (i1 + mesh.ny * i2)?
                // TODO consider changing A_exchange_field(i0+1, i1, i2) to 'local' A_exchange_field(i0, i1, i2) for x, y, z
                if(i0 == 0 && mesh.nx > 1){
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_p = a_host[util::stride(i0+1, i1, i2, mesh.nx, mesh.ny)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                        csr_ia++;
                    }
                }
                else if (i0 == mesh.nx - 1 && mesh.nx > 1){
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_m = a_host[util::stride(i0-1, i1, i2, mesh.nx, mesh.ny)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if(i0>0 && i0< mesh.nx - 1 ){
                    //i_x +u1
                    {
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_m = a_host[util::stride(i0-1, i1, i2, mesh.nx, mesh.ny)];
                    if (A_i_m != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0-1, i1, i2, im, mesh ) );
                        csr_ia++;
                    }
                    }

                    {
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_p = a_host[util::stride(i0+1, i1, i2, mesh.nx, mesh.ny)];
                    if (A_i_p != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0+1, i1, i2, im, mesh) );
                        csr_ia++;
                    }
                    }
                }

                // y
                if(i1 == 0 && mesh.ny > 1 ){
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_p = a_host[util::stride(i0, i1+1, i2, mesh.nx, mesh.ny)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if (i1 == mesh.ny - 1 && mesh.ny > 1){
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_m = a_host[util::stride(i0, i1-1, i2, mesh.nx, mesh.ny)];
                    if (A_i != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if(i1>0 && i1< mesh.ny-1){
                    {
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_m = a_host[util::stride(i0, i1-1, i2, mesh.nx, mesh.ny)];
                    if (A_i_m != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1-1, i2, im, mesh ) );
                        csr_ia++;
                    }
                    }
                    {
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_p = a_host[util::stride(i0, i1+1, i2, mesh.nx, mesh.ny)];
                    if (A_i_p != 0){
                        CSR_values.push_back( (2.* A_i)/(constants::mu0 * pow(mesh.dx, 2)) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1+1, i2, im, mesh ) );
                        csr_ia++;
                    }
                    }
                }

                // z
                if (i2 == 0 && mesh.nz > 1 ){
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_p = a_host[util::stride(i0, i1, i2+1, mesh.nx, mesh.ny)];
                    if (A_i != 0){
                        CSR_values.push_back( 2.* A_i/constants::mu0 * 1./pow(h.at(i2), 2) * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if (i2 == mesh.nz - 1 && mesh.nz > 1){
                    double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    double A_i_m = a_host[util::stride(i0, i1, i2-1, mesh.nx, mesh.ny)];
                    if (A_i != 0){
                        CSR_values.push_back( 2.* A_i/constants::mu0 * 1./pow(h.at(i2-1), 2) * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                        csr_ia++;
                    }
                }
                else if( i2 > 0 && i2 < mesh.nz - 1){
                    const double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    const double A_i_m = a_host[util::stride(i0, i1, i2-1, mesh.nx, mesh.ny)];
                    const double h_divisor = mesh.z_spacing[i2] * mesh.z_spacing[i2-1] * (1. + mesh.z_spacing[i2]/mesh.z_spacing[i2-1]);
                    if (A_i_m != 0){
                        CSR_values.push_back( 2.* A_i/constants::mu0 * (2. * h.at(i2)/h.at(i2-1))/h_divisor * 2.* A_i_m/(A_i_m + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2-1, im, mesh ) );
                        csr_ia++;
                    }
                    const double A_i_p = a_host[util::stride(i0, i1, i2+1, mesh.nx, mesh.ny)];
                    if (A_i_p != 0){
                        CSR_values.push_back( 2.* A_i/constants::mu0 * 2./h_divisor * 2.* A_i_p/(A_i_p + A_i));
                        CSR_JA.push_back( findex( i0, i1, i2+1, im, mesh ) );
                        csr_ia++;
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
    af::array result = af::sparse((dim_t) dimension, (dim_t) dimension, (dim_t) CSR_values.size(), (void*) CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if(verbose) {
        printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of CSR_matrix = %f\n", Info(), t.stop(), (float)af::sparseGetNNZ(result) / (float)result.elements());
        fflush(stdout);
    }
    return result;
}
}// namespace magnumaf
