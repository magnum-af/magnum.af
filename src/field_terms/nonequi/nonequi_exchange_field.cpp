#include "nonequi/nonequi_exchange_field.hpp"
#include "util/color_string.hpp"
#include "util/host_ptr_accessor.hpp"
#include "util/util.hpp"

namespace magnumafcpp {

af::array NonequiExchangeField::impl_H_in_Apm(const State& state) const {
    af::array exch = af::matmul(matr, af::flat(state.m));
    exch = af::moddims(exch, nemesh.nx, nemesh.ny, nemesh.nz, 3);
    if (state.Ms_field.isempty()) {
        return exch / state.Ms;
    } else {
        af::array heff = exch / state.Ms_field;
        af::replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// Get inner index (index per matrix column)
unsigned findex(unsigned i0, unsigned i1, unsigned i2, unsigned im, const NonequiMesh& mesh) {
    return i0 + mesh.nx * (i1 + mesh.ny * (i2 + mesh.nz * im));
}

af::array calc_COO_matrix(const double A_exchange, const NonequiMesh& mesh, const bool verbose) {
    printf("%s NonequiExchangeField::calc_COO_matrix unit testing not finished!\n", color_string::warning());
    fflush(stdout);
    af::timer t = af::timer::start();

    std::vector<double> h; // spacings between discretization points h = (dz[n] + dz[n+1])/2
    for (unsigned int i = 0; i < mesh.z_spacing.size() - 1; i++) {
        h.push_back((mesh.z_spacing[i] + mesh.z_spacing[i + 1]) / 2.);
    }

    const int dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_ROW;       //
    std::vector<int> COO_COL;       //
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.nz; i2++) {
            for (unsigned i1 = 0; i1 < mesh.ny; i1++) {
                for (unsigned i0 = 0; i0 < mesh.nx; i0++) {
                    const int ind = findex(i0, i1, i2, im, mesh);
                    // std::cout << ind << ", " << id << ", " << im << ", " <<
                    // i2 << ", " << i1 << ", " << i0 << std::endl; Note:
                    // skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz,
                    // 2));
                    // x
                    if (i0 == 0 && mesh.nx > 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                    }
                    if (i0 == mesh.nx - 1 && mesh.nx > 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));
                    }
                    if (i0 > 0 && i0 < mesh.nx - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));

                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                    }

                    // y
                    if (i1 == 0 && mesh.ny > 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                    }
                    if (i1 == mesh.ny - 1 && mesh.ny > 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                    }
                    if (i1 > 0 && i1 < mesh.ny - 1) {
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                    }

                    // z from ref [1]
                    if (i2 == 0 && mesh.nz > 1) { // TODO check
                        // Note: skipping f-1 term as it drops out in llg:
                        // neumann bc is assumed, which would consider fictive
                        // m[-1] with value m[0]
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(h[i2], 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                    } else if (i2 == mesh.nz - 1 && mesh.nz > 1) { // TODO check
                        // Note: skipping f+1 term as it drops out in llg:
                        // neumann bc is assumed, which would consider fictive
                        // m[n] with value m[n-1]
                        COO_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(h[i2 - 1], 2));
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2 - 1, im, mesh));
                    } else if (i2 > 0 && i2 < mesh.nz - 1) {
                        double h_divisor = h[i2] * h[i2 - 1] * (1. + h[i2] / h[i2 - 1]);
                        // f_{i-1} term
                        COO_values.push_back((2. * A_exchange / constants::mu0) * (2. * h[i2] / h[i2 - 1]) / h_divisor);
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        // f_{i+1} term
                        COO_values.push_back((2. * A_exchange / constants::mu0) * 2. / h_divisor);
                        COO_ROW.push_back(ind);
                        COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
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

// Assembly of sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array calc_COO_matrix(const af::array& A_exchange_field, const NonequiMesh& mesh, const bool verbose) {
    printf("%s NonequiExchangeField::calc_COO_matrix unit testing not finished!\n", color_string::warning());
    fflush(stdout);
    af::timer t = af::timer::start();

    std::vector<double> h; // spacings between discretization points h = (dz[n] + dz[n+1])/2
    for (unsigned int i = 0; i < mesh.z_spacing.size() - 1; i++) {
        h.push_back((mesh.z_spacing[i] + mesh.z_spacing[i + 1]) / 2.);
    }

    const int dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_ROW;       //
    std::vector<int> COO_COL;       //
    util::HostPtrAccessor<double> a_host(A_exchange_field);
    for (unsigned im = 0; im < 3; im++) {
        for (unsigned i2 = 0; i2 < mesh.nz; i2++) {
            for (unsigned i1 = 0; i1 < mesh.ny; i1++) {
                for (unsigned i0 = 0; i0 < mesh.nx; i0++) {
                    const int ind = findex(i0, i1, i2, im, mesh);
                    // Note: skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy, 2)+pow(mesh.dz,
                    // 2));
                    // x
                    if (i0 == 0 && mesh.nx > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_host[util::stride(i0 + 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                        }
                    } else if (i0 == mesh.nx - 1 && mesh.nx > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_host[util::stride(i0 - 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));
                        }
                    } else if (i0 > 0 && i0 < mesh.nx - 1) {
                        // i_x +u1
                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_m = a_host[util::stride(i0 - 1, i1, i2, mesh.nx, mesh.ny)];
                            if (A_i_m != 0) {
                                COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                     (A_i_m + A_i));
                                COO_ROW.push_back(ind);
                                COO_COL.push_back(findex(i0 - 1, i1, i2, im, mesh));
                            }
                        }

                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_p = a_host[util::stride(i0 + 1, i1, i2, mesh.nx, mesh.ny)];
                            if (A_i_p != 0) {
                                COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                     (A_i_p + A_i));
                                COO_ROW.push_back(ind);
                                COO_COL.push_back(findex(i0 + 1, i1, i2, im, mesh));
                            }
                        }
                    }

                    // y
                    if (i1 == 0 && mesh.ny > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_host[util::stride(i0, i1 + 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                        }
                    } else if (i1 == mesh.ny - 1 && mesh.ny > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_host[util::stride(i0, i1 - 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        }
                    } else if (i1 > 0 && i1 < mesh.ny - 1) {
                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_m = a_host[util::stride(i0, i1 - 1, i2, mesh.nx, mesh.ny)];
                            if (A_i_m != 0) {
                                COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_m /
                                                     (A_i_m + A_i));
                                COO_ROW.push_back(ind);
                                COO_COL.push_back(findex(i0, i1 - 1, i2, im, mesh));
                            }
                        }
                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_p = a_host[util::stride(i0, i1 + 1, i2, mesh.nx, mesh.ny)];
                            if (A_i_p != 0) {
                                COO_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_p /
                                                     (A_i_p + A_i));
                                COO_ROW.push_back(ind);
                                COO_COL.push_back(findex(i0, i1 + 1, i2, im, mesh));
                            }
                        }
                    }

                    // z
                    if (i2 == 0 && mesh.nz > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_host[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back(2. * A_i / constants::mu0 * 1. / pow(h[i2], 2) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                        }
                    } else if (i2 == mesh.nz - 1 && mesh.nz > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_host[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            COO_values.push_back(2. * A_i / constants::mu0 * 1. / pow(h[i2 - 1], 2) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        }
                    } else if (i2 > 0 && i2 < mesh.nz - 1) {
                        const double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        const double A_i_m = a_host[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];
                        const double h_divisor = h[i2] * h[i2 - 1] * (1. + h[i2] / h[i2 - 1]);
                        if (A_i_m != 0) {
                            COO_values.push_back(2. * A_i / constants::mu0 * (2. * h[i2] / h[i2 - 1]) / h_divisor * 2. *
                                                 A_i_m / (A_i_m + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        }
                        const double A_i_p = a_host[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];
                        if (A_i_p != 0) {
                            COO_values.push_back(2. * A_i / constants::mu0 * 2. / h_divisor * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(i0, i1, i2 + 1, im, mesh));
                        }
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

af::array calc_CSR_matrix(const double A_exchange, const NonequiMesh& mesh, const bool verbose) {
    printf("%s NonequiExchangeField::calc_CSR_matrix unit testing not finished!\n", color_string::warning());
    fflush(stdout);
    af::timer t = af::timer::start();

    std::vector<double> h; // spacings between discretization points h = (dz[n] + dz[n+1])/2
    for (unsigned int i = 0; i < mesh.z_spacing.size() - 1; i++) {
        h.push_back((mesh.z_spacing[i] + mesh.z_spacing[i + 1]) / 2.);
    }

    const int dimension = mesh.nx * mesh.ny * mesh.nz * 3;

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
                    const int ind = findex(i0, i1, i2, im, mesh);
                    // std::cout << ind << ", " << id << ", " << im <<
                    // ", " << i2 << ", " << i1 << ", " << i0 <<
                    // std::endl; Note: skippable due to cross product
                    // property://vmatr[findex(i0, i1, i2, im,
                    // id)]+=-6./(pow(mesh.dx, 2)+pow(mesh.dy,
                    // 2)+pow(mesh.dz, 2));
                    // x
                    if (i0 == 0 && mesh.nx > 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                        csr_ia++;
                    }
                    if (i0 == mesh.nx - 1 && mesh.nx > 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                        csr_ia++;
                    }
                    if (i0 > 0 && i0 < mesh.nx - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                        csr_ia++;

                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dx, 2));
                        CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                        csr_ia++;
                    }

                    // y
                    if (i1 == 0 && mesh.ny > 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                        csr_ia++;
                    }
                    if (i1 == mesh.ny - 1 && mesh.ny > 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        csr_ia++;
                    }
                    if (i1 > 0 && i1 < mesh.ny - 1) {
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                        csr_ia++;
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(mesh.dy, 2));
                        CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                        csr_ia++;
                    }

                    // z from ref [1]
                    if (i2 == 0 && mesh.nz > 1) { // TODO check
                        // Note: skipping f-1 term as it drops out in
                        // llg: neumann bc is assumed, which would
                        // consider fictive m[-1] with value m[0]
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(h[i2], 2));
                        CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                        csr_ia++;
                    } else if (i2 == mesh.nz - 1 && mesh.nz > 1) { // TODO check
                        // Note: skipping f+1 term as it drops out in
                        // llg: neumann bc is assumed, which would
                        // consider fictive m[n] with value m[n-1]
                        CSR_values.push_back((2. * A_exchange) / (constants::mu0)*1. / pow(h[i2 - 1], 2));
                        CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        csr_ia++;
                    } else if (i2 > 0 && i2 < mesh.nz - 1) {
                        double h_divisor = h[i2] * h[i2 - 1] * (1. + h[i2] / h[i2 - 1]);
                        // f_{i-1} term
                        CSR_values.push_back((2. * A_exchange / constants::mu0) * (2. * h[i2] / h[i2 - 1]) / h_divisor);
                        CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                        csr_ia++;
                        // f_{i+1} term
                        CSR_values.push_back((2. * A_exchange / constants::mu0) * 2. / h_divisor);
                        CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
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

// Assembly of sparse matrix for spacially varying exchange energy
// A_exchange_field
af::array calc_CSR_matrix(const af::array& A_exchange_field, const NonequiMesh& mesh, const bool verbose) {
    printf("%s NonequiExchangeField::calc_CSR_matrix unit testing not finished!\n", color_string::warning());
    fflush(stdout);
    af::timer t = af::timer::start();

    std::vector<double> h; // spacings between discretization points h = (dz[n] + dz[n+1])/2
    for (unsigned int i = 0; i < mesh.z_spacing.size() - 1; i++) {
        h.push_back((mesh.z_spacing[i] + mesh.z_spacing[i + 1]) / 2.);
    }

    const int dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> CSR_values;         // matrix values,  of length "number of elements"
    std::vector<int> CSR_IA(dimension + 1); // recursive row indices of length (n_rows + 1): IA[0] =
                                            // 0; IA[i] = IA[i-1] + (number of nonzero elements on
                                            // the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;                // comumn index of each element, hence of length
                                            // "number of elements"
    util::HostPtrAccessor<double> a_host(A_exchange_field);
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
                    // 2)+pow(mesh.dz, 2));
                    // x
                    // Note: poor indexing performace. TODO improve
                    // performance: directly accessing values with
                    // afvalue increades sp4 assembly from ~0.4 s to
                    // ~1.4 s! maybe access full host array once? is
                    // host data then in correct order for adapted
                    // findex for scalar field, i.e. i0 + mesh.nx * (i1
                    // + mesh.ny * i2)?
                    // TODO consider changing A_exchange_field(i0+1, i1,
                    // i2) to 'local' A_exchange_field(i0, i1, i2) for
                    // x, y, z
                    if (i0 == 0 && mesh.nx > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_host[util::stride(i0 + 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                            csr_ia++;
                        }
                    } else if (i0 == mesh.nx - 1 && mesh.nx > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_host[util::stride(i0 - 1, i1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                            csr_ia++;
                        }
                    } else if (i0 > 0 && i0 < mesh.nx - 1) {
                        // i_x +u1
                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_m = a_host[util::stride(i0 - 1, i1, i2, mesh.nx, mesh.ny)];
                            if (A_i_m != 0) {
                                CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_m /
                                                     (A_i_m + A_i));
                                CSR_JA.push_back(findex(i0 - 1, i1, i2, im, mesh));
                                csr_ia++;
                            }
                        }

                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_p = a_host[util::stride(i0 + 1, i1, i2, mesh.nx, mesh.ny)];
                            if (A_i_p != 0) {
                                CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dx, 2)) * 2. * A_i_p /
                                                     (A_i_p + A_i));
                                CSR_JA.push_back(findex(i0 + 1, i1, i2, im, mesh));
                                csr_ia++;
                            }
                        }
                    }

                    // y
                    if (i1 == 0 && mesh.ny > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_host[util::stride(i0, i1 + 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                            csr_ia++;
                        }
                    } else if (i1 == mesh.ny - 1 && mesh.ny > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_host[util::stride(i0, i1 - 1, i2, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                            csr_ia++;
                        }
                    } else if (i1 > 0 && i1 < mesh.ny - 1) {
                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_m = a_host[util::stride(i0, i1 - 1, i2, mesh.nx, mesh.ny)];
                            if (A_i_m != 0) {
                                CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_m /
                                                     (A_i_m + A_i));
                                CSR_JA.push_back(findex(i0, i1 - 1, i2, im, mesh));
                                csr_ia++;
                            }
                        }
                        {
                            double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                            double A_i_p = a_host[util::stride(i0, i1 + 1, i2, mesh.nx, mesh.ny)];
                            if (A_i_p != 0) {
                                CSR_values.push_back((2. * A_i) / (constants::mu0 * pow(mesh.dy, 2)) * 2. * A_i_p /
                                                     (A_i_p + A_i));
                                CSR_JA.push_back(findex(i0, i1 + 1, i2, im, mesh));
                                csr_ia++;
                            }
                        }
                    }

                    // z
                    if (i2 == 0 && mesh.nz > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_p = a_host[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back(2. * A_i / constants::mu0 * 1. / pow(h[i2], 2) * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                            csr_ia++;
                        }
                    } else if (i2 == mesh.nz - 1 && mesh.nz > 1) {
                        double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        double A_i_m = a_host[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];
                        if (A_i != 0) {
                            CSR_values.push_back(2. * A_i / constants::mu0 * 1. / pow(h[i2 - 1], 2) * 2. * A_i_m /
                                                 (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                            csr_ia++;
                        }
                    } else if (i2 > 0 && i2 < mesh.nz - 1) {
                        const double A_i = a_host[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                        const double A_i_m = a_host[util::stride(i0, i1, i2 - 1, mesh.nx, mesh.ny)];
                        const double h_divisor = h[i2] * h[i2 - 1] * (1. + h[i2] / h[i2 - 1]);
                        if (A_i_m != 0) {
                            CSR_values.push_back(2. * A_i / constants::mu0 * (2. * h[i2] / h[i2 - 1]) / h_divisor * 2. *
                                                 A_i_m / (A_i_m + A_i));
                            CSR_JA.push_back(findex(i0, i1, i2 - 1, im, mesh));
                            csr_ia++;
                        }
                        const double A_i_p = a_host[util::stride(i0, i1, i2 + 1, mesh.nx, mesh.ny)];
                        if (A_i_p != 0) {
                            CSR_values.push_back(2. * A_i / constants::mu0 * 2. / h_divisor * 2. * A_i_p /
                                                 (A_i_p + A_i));
                            CSR_JA.push_back(findex(i0, i1, i2 + 1, im, mesh));
                            csr_ia++;
                        }
                    }
                    CSR_IA[ind + 1] = CSR_IA[ind] + csr_ia;
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

NonequiExchangeField::NonequiExchangeField(NonequiMesh nemesh, double A_exchange, bool verbose, bool COO)
    : NonequiTerm(nemesh),
      matr(COO ? calc_COO_matrix(A_exchange, nemesh, verbose) : calc_CSR_matrix(A_exchange, nemesh, verbose)) {}

NonequiExchangeField::NonequiExchangeField(NonequiMesh nemesh, af::array A_exchange_field, bool verbose, bool COO)
    : NonequiTerm(nemesh), matr(COO ? calc_COO_matrix(std::move(A_exchange_field), nemesh, verbose)
                                    : calc_CSR_matrix(std::move(A_exchange_field), nemesh, verbose)) {}

// For wrapping only: constructor version taking A_exchange_field
NonequiExchangeField::NonequiExchangeField(NonequiMesh nemesh, long int A_exchange_field_ptr, bool verbose, bool COO)
    : NonequiTerm(nemesh),
      matr(COO ? calc_COO_matrix(util::pywrap::make_copy_form_py(A_exchange_field_ptr), nemesh, verbose)
               : calc_CSR_matrix(util::pywrap::make_copy_form_py(A_exchange_field_ptr), nemesh, verbose)) {}

} // namespace magnumafcpp
