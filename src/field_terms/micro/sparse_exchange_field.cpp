#include "micro/sparse_exchange_field.hpp"
#include "util/func.hpp"
#include "util/host_ptr_accessor.hpp"
#include "util/color_string.hpp"

namespace magnumafcpp {

af::array SparseExchangeField::impl_H_in_Apm(const State& state) const {
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
inline std::size_t findex(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t im, const Mesh& mesh) {
    return i0 + mesh.nx * (i1 + mesh.ny * (i2 + mesh.nz * im));
}

af::array calc_COO_matrix(const double A_exchange, const Mesh& mesh, const bool verbose) {
    af::timer t = af::timer::start();
    const unsigned dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_ROW;       //
    std::vector<int> COO_COL;       //
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.nz; i2++) {
            for (unsigned i1 = 0; i1 < mesh.ny; i1++) {
                for (unsigned i0 = 0; i0 < mesh.nx; i0++) {
                    const unsigned ind = findex(i0, i1, i2, im, mesh);
                    // std::cout << ind << ", " << id << ", " << im << ", " <<
                    // i2 << ", " << i1 << ", " << i0 << std::endl; Note:
                    // skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz,
                    // 2)); x: interaction from (i0, i1, i2) to (i0+-1, i1, i2)
                    // +x: ix, ix+1
                    if (i0 < mesh.nx - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.ny - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                    }

                    // +z: iz, iz+1
                    if (i2 < mesh.nz - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dz, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dz, 2));
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
               color_string::info(), time, time_convert,
               static_cast<double>(af::sparseGetNNZ(matr_CSR)) / static_cast<double>(matr_CSR.elements()));
        fflush(stdout);
    }
    return matr_CSR;
}

af::array calc_CSR_matrix(const double A_exchange, const Mesh& mesh, const bool verbose) {
    af::timer t = af::timer::start();
    const unsigned dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> CSR_values;         // matrix values,  of length "number of elements"
    std::vector<int> CSR_IA(dimension + 1); // recursive row indices of length (n_rows + 1): IA[0] =
                                            // 0; IA[i] = IA[i-1] + (number of nonzero elements on
                                            // the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;                // comumn index of each element, hence of length
                                            // "number of elements"
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.nz; i2++) {
            for (unsigned i1 = 0; i1 < mesh.ny; i1++) {
                for (unsigned i0 = 0; i0 < mesh.nx; i0++) {
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
                    if (i0 < mesh.nx - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                        csr_ia++;
                    }
                    // -x: ix, ix-1
                    if (i0 > 0) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                        csr_ia++;
                    }

                    // +y: iy, iy+1
                    if (i1 < mesh.ny - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                        csr_ia++;
                    }
                    // -y: iy, iy-1
                    if (i1 > 0) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        csr_ia++;
                    }

                    // +z: iz, iz+1
                    if (i2 < mesh.nz - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dz, 2));
                        CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                        csr_ia++;
                    }
                    // -z: iz, iz-1
                    if (i2 > 0) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0) / pow(mesh.dz, 2));
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
               color_string::info(), t.stop(),
               static_cast<double>(af::sparseGetNNZ(result)) / static_cast<double>(result.elements()));
    return result;
}

inline auto lapA(const double A_i, const double A_pm, const double dxyz) {
    return 2. * A_i / (constants::mu0 * pow(dxyz, 2)) * 2. * A_pm / (A_pm + A_i);
}

// Assembly of sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array calc_CSR_matrix(const af::array& A_exchange_field, const Mesh& mesh, const bool verbose) {
    printf("%s SparseExchangeField::calc_CSR_matrix unit testing not finished!\n", color_string::warning());
    fflush(stdout);
    af::timer t = af::timer::start();
    const unsigned dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> CSR_values;         // matrix values,  of length "number of elements"
    std::vector<int> CSR_IA(dimension + 1); // recursive row indices of length (n_rows + 1): IA[0] =
                                            // 0; IA[i] = IA[i-1] + (number of nonzero elements on
                                            // the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;                // comumn index of each element, hence of length
                                            // "number of elements"
    util::HostPtrAccessor<double> a_host(A_exchange_field);
    for (std::size_t im = 0; im < 3; im++) {
        for (std::size_t i2 = 0; i2 < mesh.nz; i2++) {
            for (std::size_t i1 = 0; i1 < mesh.ny; i1++) {
                for (std::size_t i0 = 0; i0 < mesh.nx; i0++) {
                    int csr_ia = 0; // counter for SCR_IA
                    const double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];

                    if (A_i != 0.0) {
                        double m_self{0.0};

                        auto push_CSR_val = [&a_host, &CSR_values, &CSR_JA, &csr_ia, &m_self, mesh, A_i,
                                             im](const auto ix, const auto iy, const auto iz, const auto dxyz) {
                            const double A_i_pm = a_host[util::stride(ix, iy, iz, mesh.nx, mesh.ny)];
                            const auto lapA_val = lapA(A_i, A_i_pm, dxyz);
                            CSR_values.push_back(lapA_val);
                            CSR_JA.push_back(findex(ix, iy, iz, im, mesh));
                            csr_ia++;
                            m_self -= lapA_val;
                        };

                        // TODO set middle element of laplace?
                        // std::cout << ind << ", " << id << ", " << im <<
                        // ", " << i2 << ", " << i1 << ", " << i0 <<
                        // std::endl; Note: skippable due to cross product
                        // property://vmatr[findex(i0, i1, i2, im,
                        // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy,
                        // 2)+pow(mesh.dz, 2));

                        // +x: ix, ix+1
                        if (i0 < mesh.nx - 1) {
                            push_CSR_val(i0 + 1, i1, i2, mesh.dx);
                        }
                        // TODO: neuman BC
                        // else {
                        //     push_CSR_val(i0, i1, i2, mesh.dx);
                        // }
                        // -x: ix, ix-1
                        if (i0 > 0) {
                            push_CSR_val(i0 - 1, i1, i2, mesh.dx);
                        }

                        // +y: iy, iy+1
                        if (i1 < mesh.ny - 1) {
                            push_CSR_val(i0, i1 + 1, i2, mesh.dy);
                        }
                        // -y: iy, iy-1
                        if (i1 > 0) {
                            push_CSR_val(i0, i1 - 1, i2, mesh.dy);
                        }

                        // +z: iz, iz+1
                        if (i2 < mesh.nz - 1) {
                            push_CSR_val(i0, i1, i2 + 1, mesh.dz);
                        }
                        // -z: iz, iz-1
                        if (i2 > 0) {
                            push_CSR_val(i0, i1, i2 - 1, mesh.dz);
                        }

                        // adding m_i
                        CSR_values.push_back(m_self);
                        CSR_JA.push_back(findex(i0, i1, i2, im, mesh));
                        csr_ia++;
                    }

                    const auto ind = findex(i0, i1, i2, im, mesh);
                    CSR_IA[ind + 1] = CSR_IA[ind] + csr_ia; // Called at each iteration, don't move into if
                }
            }
        }
    }
    af::array result = af::sparse((dim_t)dimension, (dim_t)dimension, (dim_t)CSR_values.size(),
                                  (void*)CSR_values.data(), CSR_IA.data(), CSR_JA.data(), f64);
    if (verbose) {
        printf("%s Initialized sparse exchange matrix in %f [s]. Sparsity of "
               "CSR_matrix = %f\n",
               color_string::info(), t.stop(),
               static_cast<double>(af::sparseGetNNZ(result)) / static_cast<double>(result.elements()));
        fflush(stdout);
    }
    return result;
}

// Assembly of COO sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array calc_COO_matrix(const af::array& A_exchange_field, const Mesh& mesh, const bool verbose) {
    printf("%s SparseExchangeField::calc_COO_matrix unit testing not finished!\n", color_string::warning());
    fflush(stdout);
    af::timer t = af::timer::start();
    const std::size_t dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_COL;       // column indices
    std::vector<int> COO_ROW;       // row indices
    util::HostPtrAccessor<double> a_raw(A_exchange_field);
    for (std::size_t im = 0; im < 3; im++) {
        for (std::size_t i2 = 0; i2 < mesh.nz; i2++) {
            for (std::size_t i1 = 0; i1 < mesh.ny; i1++) {
                for (std::size_t i0 = 0; i0 < mesh.nx; i0++) {
                    const std::size_t ind = findex(i0, i1, i2, im, mesh);
                    const double A_i = a_raw[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    if (A_i != 0.0) {
                        // m_i element, will be added up in lambda to -6*...
                        double m_self{0.0};

                        // TODO set middle element of laplace?
                        // std::cout << ind << ", " << im << ", " << i2 << ", " << i1 << ", " << i0 << std::endl;
                        // Note:
                        // skippable due to cross product
                        // property://vmatr[findex(i0, i1, i2, im,
                        // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz,
                        // 2));

                        auto push_val = [im, &a_raw, &COO_values, &COO_ROW, &COO_COL, mesh, ind, A_i,
                                         &m_self](const auto ix, const auto iy, const auto iz, const auto dxyz) {
                            const double A_i_pm = a_raw[util::stride(ix, iy, iz, mesh.nx, mesh.ny)];
                            const auto lapA_val = lapA(A_i, A_i_pm, dxyz);
                            COO_values.push_back(lapA_val);
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(ix, iy, iz, im, mesh));
                            m_self -= lapA_val;
                        };

                        // +x: ix, ix+1
                        if (i0 < mesh.nx - 1) {
                            push_val(i0 + 1, i1, i2, mesh.dx);
                        }
                        // -x: ix, ix-1
                        if (i0 > 0) {
                            push_val(i0 - 1, i1, i2, mesh.dx);
                        }

                        // +y: iy, iy+1
                        if (i1 < mesh.ny - 1) {
                            push_val(i0, i1 + 1, i2, mesh.dy);
                        }
                        // -y: iy, iy-1
                        if (i1 > 0) {
                            push_val(i0, i1 - 1, i2, mesh.dy);
                        }

                        // +z: iz, iz+1
                        if (i2 < mesh.nz - 1) {
                            push_val(i0, i1, i2 + 1, mesh.dz);
                        }
                        // -z: iz, iz-1
                        if (i2 > 0) {
                            push_val(i0, i1, i2 - 1, mesh.dz);
                        }

                        // adding m_i
                        COO_values.push_back(m_self);
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2, im, mesh));
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
               color_string::info(), time, time_convert,
               static_cast<double>(af::sparseGetNNZ(matr_CSR)) / static_cast<double>(matr_CSR.elements()));
        fflush(stdout);
    }
    return matr_CSR;
}

SparseExchangeField::SparseExchangeField(double A_exchange, Mesh mesh, bool verbose, bool COO)
    : matr(COO ? calc_COO_matrix(A_exchange, mesh, verbose) : calc_CSR_matrix(A_exchange, mesh, verbose)) {}

SparseExchangeField::SparseExchangeField(const af::array& A_exchange_field, Mesh mesh, bool verbose, bool COO)
    : matr(COO ? calc_COO_matrix(A_exchange_field, mesh, verbose) : calc_CSR_matrix(A_exchange_field, mesh, verbose)) {}

// For wrapping only: constructor version taking A_exchange_field
SparseExchangeField::SparseExchangeField(long int A_exchange_field_ptr, Mesh mesh, bool verbose)
    : matr(calc_CSR_matrix(*(new af::array(*((void**)A_exchange_field_ptr))), mesh, verbose)) {}

} // namespace magnumafcpp
