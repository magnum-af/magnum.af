#include "micro_exch_sparse.hpp"
#include "util/func.hpp"
#include "util/misc.hpp"

namespace magnumafcpp {

SparseExchangeField::SparseExchangeField(double A_exchange, Mesh mesh, bool verbose, bool COO)
    : matr(COO ? calc_COO_matrix(A_exchange, mesh, verbose) : calc_CSR_matrix(A_exchange, mesh, verbose)) {}

SparseExchangeField::SparseExchangeField(const af::array& A_exchange_field, Mesh mesh, bool verbose, bool COO)
    : matr(COO ? calc_COO_matrix(A_exchange_field, mesh, verbose) : calc_CSR_matrix(A_exchange_field, mesh, verbose)) {}

// For wrapping only: constructor version taking A_exchange_field
SparseExchangeField::SparseExchangeField(long int A_exchange_field_ptr, Mesh mesh, bool verbose)
    : matr(calc_CSR_matrix(*(new af::array(*((void**)A_exchange_field_ptr))), mesh, verbose)) {}

af::array SparseExchangeField::h(const State& state) {
    af::timer aftimer = af::timer::start();
    af::array exch = af::matmul(matr, af::flat(state.m));
    exch = af::moddims(exch, state.mesh.n0, state.mesh.n1, state.mesh.n2, 3);
    if (state.afsync)
        af::sync();
    af_time += af::timer::stop(aftimer);
    if (state.Ms_field.isempty()) {
        return exch / state.Ms;
    } else {
        af::array heff = exch / state.Ms_field;
        replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// Get inner index (index per matrix column)
unsigned SparseExchangeField::findex(unsigned i0, unsigned i1, unsigned i2, unsigned im, const Mesh& mesh) {
    return i0 + mesh.n0 * (i1 + mesh.n1 * (i2 + mesh.n2 * im));
}

af::array SparseExchangeField::calc_COO_matrix(const double A_exchange, const Mesh& mesh, const bool verbose) {
    af::timer t;
    if (verbose)
        af::timer::start();
    const unsigned dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;

    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_ROW;       //
    std::vector<int> COO_COL;       //
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.n2; i2++) {
            for (unsigned i1 = 0; i1 < mesh.n1; i1++) {
                for (unsigned i0 = 0; i0 < mesh.n0; i0++) {
                    const unsigned ind = findex(i0, i1, i2, im, mesh);
                    // std::cout << ind << ", " << id << ", " << im << ", " <<
                    // i2 << ", " << i1 << ", " << i0 << std::endl; Note:
                    // skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz,
                    // 2)); x: interaction from (i0, i1, i2) to (i0+-1, i1, i2)
                    // +x: ix, ix+1
                    if (i0 < mesh.n0 - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.n1 - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                    }

                    // +z: iz, iz+1
                    if (i2 < mesh.n2 - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dz, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dz, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2 - 1, im, mesh));
                    }
                }
            }
        }
    }
    af::array matr_COO = af::sparse((dim_t)dimension, (dim_t)dimension, af::array(COO_values.size(), COO_values.data()),
                                    af::array(COO_ROW.size(), COO_ROW.data()),
                                    af::array(COO_COL.size(), COO_COL.data()), AF_STORAGE_COO);
    double time = t.stop();

    af::timer timer_convert = af::timer::start();
    af::array matr_CSR = af::sparseConvertTo(matr_COO, AF_STORAGE_CSR);
    double time_convert = timer_convert.stop();

    if (verbose) {
        printf("%s Initialized sparse COO exchange matrix in %f [s]. Converted "
               "COO to CSR in %f [s]. Sparsity = %f\n",
               Info(), time, time_convert,
               static_cast<double>(af::sparseGetNNZ(matr_CSR)) / static_cast<double>(matr_CSR.elements()));
        fflush(stdout);
    }
    return matr_CSR;
}

af::array SparseExchangeField::calc_CSR_matrix(const double A_exchange, const Mesh& mesh, const bool verbose) {
    af::timer t;
    if (verbose)
        af::timer::start();
    const unsigned dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;

    std::vector<double> CSR_values;         // matrix values,  of length "number of elements"
    std::vector<int> CSR_IA(dimension + 1); // recursive row indices of length (n_rows + 1): IA[0] =
                                            // 0; IA[i] = IA[i-1] + (number of nonzero elements on
                                            // the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;                // comumn index of each element, hence of length
                                            // "number of elements"
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.n2; i2++) {
            for (unsigned i1 = 0; i1 < mesh.n1; i1++) {
                for (unsigned i0 = 0; i0 < mesh.n0; i0++) {
                    int csr_ia = 0; // counter for SCR_IA
                    const unsigned ind = findex(i0, i1, i2, im, mesh);
                    // std::cout << ind << ", " << id << ", " << im <<
                    // ", " << i2 << ", " << i1 << ", " << i0 <<
                    // std::endl; Note: skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy,
                    // 2)+pow(mesh.dz, 2)); x: interaction from (i0, i1,
                    // i2) to (i0+-1, i1, i2)
                    // +x: ix, ix+1
                    if (i0 < mesh.n0 - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                        csr_ia++;
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                        csr_ia++;
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.n1 - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                        csr_ia++;
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        csr_ia++;
                    }

                    // +z: iz, iz+1
                    if (i2 < mesh.n2 - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dz, 2));
                        CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                        csr_ia++;
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dz, 2));
                        CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        csr_ia++;
                    }
                    CSR_IA[ind + 1] = CSR_IA[ind] + csr_ia;
                }
            }
        }
    }

    af::array result = af::sparse((dim_t)dimension, (dim_t)dimension, (dim_t)CSR_values.size(),
                                  (void*)CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if (verbose)
        printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of "
               "CSR_matrix = %f\n",
               Info(), t.stop(), static_cast<double>(af::sparseGetNNZ(matr)) / static_cast<double>(matr.elements()));
    return result;
}

// Assembly of sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array SparseExchangeField::calc_CSR_matrix(const af::array& A_exchange_field, const Mesh& mesh,
                                               const bool verbose) {
    printf("%s SparseExchangeField::calc_CSR_matrix unit testing not finished!\n", Warning());
    fflush(stdout);
    af::timer t;
    if (verbose)
        af::timer::start();
    const unsigned dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;

    std::vector<double> CSR_values;         // matrix values,  of length "number of elements"
    std::vector<int> CSR_IA(dimension + 1); // recursive row indices of length (n_rows + 1): IA[0] =
                                            // 0; IA[i] = IA[i-1] + (number of nonzero elements on
                                            // the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;                // comumn index of each element, hence of length
                                            // "number of elements"
    double* a_host = NULL;
    a_host = A_exchange_field.host<double>();
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.n2; i2++) {
            for (unsigned i1 = 0; i1 < mesh.n1; i1++) {
                for (unsigned i0 = 0; i0 < mesh.n0; i0++) {
                    int csr_ia = 0; // counter for SCR_IA
                    const unsigned ind = findex(i0, i1, i2, im, mesh);
                    // std::cout << ind << ", " << id << ", " << im <<
                    // ", " << i2 << ", " << i1 << ", " << i0 <<
                    // std::endl; Note: skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy,
                    // 2)+pow(mesh.dz, 2));
                    // +x: ix, ix+1
                    if (i0 < mesh.n0 - 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_host[util::stride(i0 + 1, i1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                            csr_ia++;
                        }
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_host[util::stride(i0 - 1, i1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                            csr_ia++;
                        }
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.n1 - 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_host[util::stride(i0, i1 + 1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                            csr_ia++;
                        }
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_host[util::stride(i0, i1 - 1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                            csr_ia++;
                        }
                    }

                    // +z: iz, iz+1
                    if (i2 < mesh.n2 - 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_host[util::stride(i0, i1, i2 + 1, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dz, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                            csr_ia++;
                        }
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_host[util::stride(i0, i1, i2 - 1, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dz, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                            csr_ia++;
                        }
                    }
                    CSR_IA[ind + 1] = CSR_IA[ind] + csr_ia;
                }
            }
        }
    }
    af::freeHost(a_host);
    af::array result = af::sparse((dim_t)dimension, (dim_t)dimension, (dim_t)CSR_values.size(),
                                  (void*)CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if (verbose) {
        printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of "
               "CSR_matrix = %f\n",
               Info(), t.stop(),
               static_cast<double>(af::sparseGetNNZ(result)) / static_cast<double>(result.elements()));
        fflush(stdout);
    }
    return result;
}

// Assembly of COO sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array SparseExchangeField::calc_COO_matrix(const af::array& A_exchange_field, const Mesh& mesh,
                                               const bool verbose) {
    printf("%s SparseExchangeField::calc_COO_matrix unit testing not finished!\n", Warning());
    fflush(stdout);
    af::timer t;
    if (verbose)
        af::timer::start();
    const unsigned dimension = mesh.n0 * mesh.n1 * mesh.n2 * 3;

    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_COL;       // column indices
    std::vector<int> COO_ROW;       // row indices
    double* a_raw = NULL;
    a_raw = A_exchange_field.host<double>();
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.n2; i2++) {
            for (unsigned i1 = 0; i1 < mesh.n1; i1++) {
                for (unsigned i0 = 0; i0 < mesh.n0; i0++) {
                    const int ind = findex(i0, i1, i2, im, mesh);
                    // std::cout << ind << ", " << id << ", " << im << ", " <<
                    // i2 << ", " << i1 << ", " << i0 << std::endl; Note:
                    // skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz,
                    // 2));
                    // +x: ix, ix+1
                    if (i0 < mesh.n0 - 1) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_raw[util::stride(i0 + 1, i1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                        }
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_raw[util::stride(i0 - 1, i1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));
                        }
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.n1 - 1) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_raw[util::stride(i0, i1 + 1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                        }
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_raw[util::stride(i0, i1 - 1, i2, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        }
                    }

                    // +z: iz, iz+1
                    if (i2 < mesh.n2 - 1) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_p = a_raw[util::stride(i0, i1, i2 + 1, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dz, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                        }
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.n0, mesh.n1)];
                        double A_i_m = a_raw[util::stride(i0, i1, i2 - 1, mesh.n0, mesh.n1)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dz, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        }
                    }
                }
            }
        }
    }
    af::freeHost(a_raw);
    af::array matr_COO = af::sparse((dim_t)dimension, (dim_t)dimension, af::array(COO_values.size(), COO_values.data()),
                                    af::array(COO_ROW.size(), COO_ROW.data()),
                                    af::array(COO_COL.size(), COO_COL.data()), AF_STORAGE_COO);
    double time = t.stop();

    af::timer timer_convert = af::timer::start();
    af::array matr_CSR = af::sparseConvertTo(matr_COO, AF_STORAGE_CSR);
    double time_convert = timer_convert.stop();

    if (verbose) {
        printf("%s Initialized sparse COO exchange matrix in %f [s]. Converted "
               "COO to CSR in %f [s]. Sparsity = %f\n",
               Info(), time, time_convert,
               static_cast<double>(af::sparseGetNNZ(matr_CSR)) / static_cast<double>(matr_CSR.elements()));
        fflush(stdout);
    }
    return matr_CSR;
}
} // namespace magnumafcpp
