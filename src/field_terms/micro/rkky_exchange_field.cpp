#include "micro/rkky_exchange_field.hpp"
#include "util/color_string.hpp"
#include "util/func.hpp"
#include "util/host_ptr_accessor.hpp"
#include "util/util.hpp"

namespace magnumafcpp {

RKKYExchangeField::RKKYExchangeField(RKKY_values rkky_values, Exchange_values exchange_values, Mesh mesh,
                                     const af::array& rkky_indices, bool verbose, bool COO)
    : matr(COO ? calc_COO_matrix(rkky_values.get(), exchange_values.get(), mesh, rkky_indices, verbose)
               : calc_CSR_matrix(rkky_values.get(), exchange_values.get(), mesh, rkky_indices, verbose)) {}

RKKYExchangeField::RKKYExchangeField(long int rkky_values, long int exchange_values, Mesh mesh, long int rkky_indices,
                                     bool verbose)
    : matr(calc_COO_matrix(util::pywrap::make_copy_form_py(rkky_values), util::pywrap::make_copy_form_py(exchange_values), mesh,
                           util::pywrap::make_copy_form_py(rkky_indices), verbose)) {}

af::array RKKYExchangeField::impl_H_in_Apm(const State& state) const {
    af::array exch = af::matmul(matr, af::flat(state.m));
    exch = af::moddims(exch, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3);
    if (state.Ms_field.isempty()) {
        return exch / state.Ms;
    } else {
        af::array heff = exch / state.Ms_field;
        af::replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// Get inner index (index per matrix column)
int RKKYExchangeField::findex(unsigned i0, unsigned i1, unsigned i2, unsigned im, const Mesh& mesh) const {
    return static_cast<int>(i0 + mesh.nx * (i1 + mesh.ny * (i2 + mesh.nz * im)));
}

// Assembly of sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array RKKYExchangeField::calc_CSR_matrix(const af::array& RKKY_field, const af::array& A_exchange_field,
                                             const Mesh& mesh, const af::array& rkky_indices,
                                             const bool verbose) const {
    printf("%s RKKYExchangeField::calc_CSR_matrix unit testing not finished!\n", color_string::warning());
    fflush(stdout);
    af::timer t = af::timer::start();
    const unsigned dimension = mesh.nx * mesh.ny * mesh.nz * 3;
    // matrix values,  of length "number of elements"
    std::vector<double> CSR_values;
    std::vector<int> CSR_IA(dimension + 1); // recursive row indices of length (n_rows + 1): IA[0] =
                                            // 0; IA[i] = IA[i-1] + (number of nonzero elements on
                                            // the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;                // comumn index of each element, hence of length
                                            // "number of elements"
    util::HostPtrAccessor<double> a_raw(A_exchange_field);
    util::HostPtrAccessor<double> rkky_raw(RKKY_field);
    util::HostPtrAccessor<unsigned> rkky_indices_raw(rkky_indices); // .ptr is set to nullptr if a.isempty()

    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.nz; i2++) {
            for (unsigned i1 = 0; i1 < mesh.ny; i1++) {
                for (unsigned i0 = 0; i0 < mesh.nx; i0++) {
                    unsigned csr_ia = 0; // counter for SCR_IA
                    const int ind = findex(i0, i1, i2, im, mesh);
                    // Note: skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy,
                    // 2)+pow(mesh.dz, 2));
                    // +x: ix, ix+1
                    if (i0 < mesh.nx - 1) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_raw[util::stride(i0 + 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                            csr_ia++;
                        }
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_raw[util::stride(i0 - 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                            csr_ia++;
                        }
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.ny - 1) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_raw[util::stride(i0, i1 + 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                            csr_ia++;
                        }
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_raw[util::stride(i0, i1 - 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                            csr_ia++;
                        }
                    }

                    // +z: iz, iz+1
                    if (i2 < mesh.nz - 1) {
                        double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double RKKY_i_p = rkky_raw[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];

                        const unsigned int RKKY_index_i =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)] : 0;
                        const unsigned int RKKY_index_i_p =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)] : 0;
                        // std::cout << "rkkyindex = " << RKKY_index_i
                        // << "and" << RKKY_index_i_p << std::endl;
                        // Preferring RKKY over exch vals
                        if ((RKKY_index_i == RKKY_index_i_p) && (RKKY_i != 0) && (RKKY_i_p != 0)) {
                            // assuming rkky jump condition equal to
                            // exch jump
                            CSR_values.push_back((2. * RKKY_i) / (constants::mu0 * pow(mesh.dz, 2)));
                            CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                            csr_ia++;
                        } else {
                            double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_p = a_raw[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];
                            if (A_i != 0) {
                                CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dz, 2)) * 2. * A_i_p /
                                                     (A_i_p + A_i));
                                CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                                csr_ia++;
                            }
                        }
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double RKKY_i_m = rkky_raw[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];

                        const unsigned int RKKY_index_i =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)] : 0;
                        const unsigned int RKKY_index_i_m =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)] : 0;

                        if ((RKKY_index_i == RKKY_index_i_m) && (RKKY_i != 0) && (RKKY_i_m != 0)) {
                            // assuming rkky jump condition equal to
                            // exch jump
                            CSR_values.push_back((2. * RKKY_i) / (constants::mu0 * pow(mesh.dz, 2)));
                            CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                            csr_ia++;
                        } else {
                            double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_m = a_raw[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];
                            if (A_i != 0) {
                                CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dz, 2)) * 2. * A_i_m /
                                                     (A_i_m + A_i));
                                CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                                csr_ia++;
                            }
                        }
                    }
                    CSR_IA[ind + 1] = CSR_IA[ind] + csr_ia;
                }
            }
        }
    }

    // for (auto const& value: CSR_IA){
    //    std::cout << "CSR_IA=" << value << std::endl;
    //}
    // for (auto const& value: CSR_JA){
    //    std::cout << "CSR_JA=" << value << std::endl;
    //}
    af::array result = af::sparse((dim_t)dimension, (dim_t)dimension, (dim_t)CSR_values.size(),
                                  (void*)CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if (verbose) {
        printf("%s Initialized sparse CSR RKKY-exchange matrix in %f [s]. "
               "Sparsity of CSR_matrix = %f\n",
               color_string::info(), t.stop(),
               static_cast<double>(af::sparseGetNNZ(result)) / static_cast<double>(result.elements()));
        fflush(stdout);
    }
    return result;
}

// Assembly of sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array RKKYExchangeField::calc_COO_matrix(const af::array& RKKY_field, const af::array& A_exchange_field,
                                             const Mesh& mesh, const af::array& rkky_indices,
                                             const bool verbose) const {
    std::cout << color_string::info() << " Starting RKKYExchangeField sparse matrix setup" << std::endl;
    af::timer t = af::timer::start();
    const unsigned dimension = mesh.nx * mesh.ny * mesh.nz * 3;
    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_COL;
    std::vector<int> COO_ROW;

    util::HostPtrAccessor<double> a_raw(A_exchange_field);
    util::HostPtrAccessor<double> rkky_raw(RKKY_field);
    util::HostPtrAccessor<unsigned> rkky_indices_raw(rkky_indices); // .ptr is set to nullptr if a.isempty()
    // NOTE aborts program//#pragma omp parallel for
    // consider removing im loop, but tiling with sparse is not supported.
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.nz; i2++) {
            for (unsigned i1 = 0; i1 < mesh.ny; i1++) {
                for (unsigned i0 = 0; i0 < mesh.nx; i0++) {
                    const int ind = findex(i0, i1, i2, im, mesh);
                    // +x: ix, ix+1
                    if (i0 < mesh.nx - 1) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_raw[util::stride(i0 + 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                        }
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_raw[util::stride(i0 - 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));
                        }
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.ny - 1) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_raw[util::stride(i0, i1 + 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                        }
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_raw[util::stride(i0, i1 - 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        }
                    }

                    // z
                    // +z: iz, iz+1
                    if (i2 < mesh.nz - 1) {
                        double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double RKKY_i_p = rkky_raw[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];

                        const unsigned int RKKY_index_i =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)] : 0;
                        const unsigned int RKKY_index_i_p =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)] : 0;
                        // std::cout << "rkkyindex = " << RKKY_index_i << "and"
                        // << RKKY_index_i_p << std::endl;
                        // Preferring RKKY over exch vals:
                        if ((RKKY_index_i == RKKY_index_i_p) && (RKKY_i != 0) && (RKKY_i_p != 0)) {
                            COO_values.push_back((2. * RKKY_i) / (constants::mu0 * pow(mesh.dz, 2)));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                        } else {
                            double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_p = a_raw[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];
                            if (A_i != 0) {
                                COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dz, 2)) * 2. * A_i_p /
                                                     (A_i_p + A_i));
                                COO_ROW.push_back(ind);
                                COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                            }
                        }
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        double RKKY_i = rkky_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double RKKY_i_m = rkky_raw[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];

                        const unsigned int RKKY_index_i =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)] : 0;
                        const unsigned int RKKY_index_i_m =
                            rkky_indices_raw ? rkky_indices_raw[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)] : 0;

                        if ((RKKY_index_i == RKKY_index_i_m) && (RKKY_i != 0) && (RKKY_i_m != 0)) {
                            COO_values.push_back((2. * RKKY_i) / (constants::mu0 * pow(mesh.dz, 2)));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        } else {
                            double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_m = a_raw[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];
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
    }
    af::array matr_COO = af::sparse((dim_t)dimension, (dim_t)dimension,
                                    af::array(static_cast<dim_t>(COO_values.size()), COO_values.data()),
                                    af::array(static_cast<dim_t>(COO_ROW.size()), COO_ROW.data()),
                                    af::array(static_cast<dim_t>(COO_COL.size()), COO_COL.data()), AF_STORAGE_COO);
    double time = t.stop();

    af::timer timer_convert = af::timer::start();
    af::array matr_CSR = af::sparseConvertTo(matr_COO, AF_STORAGE_CSR);
    double time_convert = timer_convert.stop();

    if (verbose) {
        printf("%s Finished sparse matrix setup in %f [s]. Initialized COO matrix "
               "in %f [s]. Converted COO to CSR format in %f [s]. Sparsity = %f\n",
               color_string::info(), time + time_convert, time, time_convert,
               static_cast<double>(af::sparseGetNNZ(matr_CSR)) / static_cast<double>(matr_CSR.elements()));
        fflush(stdout);
    }

    return matr_CSR;
}

} // namespace magnumafcpp
