#include "micro/exchange_field_pbc.hpp"
#include "util/func.hpp"
#include "util/host_ptr_accessor.hpp"
#include "util/misc.hpp"
#include <optional>

namespace magnumafcpp {

af::array ExchangeFieldPBC::impl_H_in_Apm(const State& state) const {
    af::array exch = af::matmul(matr, af::flat(state.m));
    exch = af::moddims(exch, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3);
    if (state.Ms_field.isempty()) {
        return exch / state.Ms;
    } else {
        af::array heff = exch / state.Ms_field;
        replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// Get inner index (index per matrix column)
inline std::size_t findex(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t im, const Mesh& mesh) {
    return i0 + mesh.nx * (i1 + mesh.ny * (i2 + mesh.nz * im));
}

inline auto lapA(const double A_i, const double A_pm, const double dxyz) {
    return 2. * A_i / (constants::mu0 * pow(dxyz, 2)) * 2. * A_pm / (A_pm + A_i);
}

// Assembly of sparse matrix for spacially varying exchange energy
// A_exchange_field
// af::array PBC_CSR_matrix(const af::array& A_exchange_field, const Mesh& mesh, const bool verbose) {
af::array PBC_CSR_matrix(const std::optional<af::array>& A_exchange_field, const std::optional<double> A_scalar,
                         const Mesh& mesh, const bool verbose) {
    // printf("%s ExchangeFieldPBC::PBC_CSR_matrix unit testing not finished!\n", Warning());
    // fflush(stdout);
    af::timer t = af::timer::start();
    const unsigned dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> CSR_values;         // matrix values,  of length "number of elements"
    std::vector<int> CSR_IA(dimension + 1); // recursive row indices of length (n_rows + 1): IA[0] =
                                            // 0; IA[i] = IA[i-1] + (number of nonzero elements on
                                            // the i-1-th row in the original matrix)
    std::vector<int> CSR_JA;                // comumn index of each element, hence of length
                                            // "number of elements"
    std::optional<util::HostPtrAccessor<double>> a_raw;
    if (A_exchange_field) {
        a_raw = A_exchange_field;
    }
    for (std::size_t im = 0; im < 3; im++) {
        for (std::size_t i2 = 0; i2 < mesh.nz; i2++) {
            for (std::size_t i1 = 0; i1 < mesh.ny; i1++) {
                for (std::size_t i0 = 0; i0 < mesh.nx; i0++) {
                    int csr_ia = 0; // counter for SCR_IA
                    const auto ind = findex(i0, i1, i2, im, mesh);
                    const double A_i =
                        A_scalar ? A_scalar.value() : a_raw.value()[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    if (A_i != 0.0) {
                        double m_self{0.0}; // m_i element, will be added up in lambda to -6*...
                        auto push_CSR_val = [&m_self, im, &a_raw, &CSR_values, &CSR_JA, &csr_ia, mesh, A_i, A_scalar](
                                                const auto ix_in, const auto iy_in, const auto iz_in, const auto dxyz) {
                            // std::cout << ix_in << " " << (ix_in + mesh.nx) % mesh.nx << std::endl;
                            const auto ix = (ix_in + mesh.nx) % mesh.nx;
                            const auto iy = (iy_in + mesh.ny) % mesh.ny;
                            const auto iz = (iz_in + mesh.nz) % mesh.nz;
                            const double A_i_pm =
                                A_scalar ? A_scalar.value() : a_raw.value()[util::stride(ix, iy, iz, mesh.nx, mesh.ny)];
                            const auto lapA_val = lapA(A_i, A_i_pm, dxyz);
                            CSR_values.push_back(lapA_val);
                            CSR_JA.push_back(findex(ix, iy, iz, im, mesh));
                            csr_ia++;
                            m_self -= lapA_val;
                        };

                        push_CSR_val(i0 + 1, i1, i2, mesh.dx); // +x: ix, ix+1
                        push_CSR_val(i0 - 1, i1, i2, mesh.dx); // -x: ix, ix-1
                        push_CSR_val(i0, i1 + 1, i2, mesh.dy); // +y: iy, iy+1
                        push_CSR_val(i0, i1 - 1, i2, mesh.dy); // -y: iy, iy-1
                        push_CSR_val(i0, i1, i2 + 1, mesh.dz); // +z: iz, iz+1
                        push_CSR_val(i0, i1, i2 - 1, mesh.dz); // -z: iz, iz-1

                        // add m_i element
                        CSR_values.push_back(m_self);
                        CSR_JA.push_back(ind);
                        csr_ia++;
                    }

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
               Info(), t.stop(),
               static_cast<double>(af::sparseGetNNZ(result)) / static_cast<double>(result.elements()));
        fflush(stdout);
    }
    return result;
}

// Assembly of COO sparse matrix for spacially varying exchange energy
// A_exchange_field
// Expects either first or second optional, not both
af::array PBC_COO_matrix(const std::optional<af::array>& A_exchange_field, const std::optional<double> A_scalar,
                         const Mesh& mesh, const bool verbose) {
    // af::array PBC_COO_matrix(const af::array& A_exchange_field, const Mesh& mesh, const bool verbose) {
    // printf("%s ExchangeFieldPBC::PBC_COO_matrix unit testing not finished!\n", Warning());
    // fflush(stdout);
    af::timer t = af::timer::start();
    const std::size_t dimension = mesh.nx * mesh.ny * mesh.nz * 3;

    std::vector<double> COO_values; // matrix values,  of length "number of elements"
    std::vector<int> COO_COL;       // column indices
    std::vector<int> COO_ROW;       // row indices

    std::optional<util::HostPtrAccessor<double>> a_raw;
    if (A_exchange_field) {
        a_raw = A_exchange_field;
    }

    for (std::size_t im = 0; im < 3; im++) {
        for (std::size_t i2 = 0; i2 < mesh.nz; i2++) {
            for (std::size_t i1 = 0; i1 < mesh.ny; i1++) {
                for (std::size_t i0 = 0; i0 < mesh.nx; i0++) {
                    const std::size_t ind = findex(i0, i1, i2, im, mesh);
                    const double A_i =
                        A_scalar ? A_scalar.value() : a_raw.value()[util::stride(i0, i1, i2, mesh.nx, mesh.ny)];
                    if (A_i != 0.0) {
                        double m_self{0.0}; // m_i element, will be added up in lambda to -6*...
                        auto push_val = [&m_self, im, &a_raw, &COO_values, &COO_ROW, &COO_COL, mesh, ind, A_i,
                                         A_scalar](const auto ix_in, const auto iy_in, const auto iz_in,
                                                   const auto dxyz) {
                            // std::cout << ix_in << " " << (ix_in + mesh.nx) % mesh.nx << std::endl;
                            const auto ix = (ix_in + mesh.nx) % mesh.nx;
                            const auto iy = (iy_in + mesh.ny) % mesh.ny;
                            const auto iz = (iz_in + mesh.nz) % mesh.nz;
                            const double A_i_pm =
                                A_scalar ? A_scalar.value() : a_raw.value()[util::stride(ix, iy, iz, mesh.nx, mesh.ny)];
                            const auto lapA_val = lapA(A_i, A_i_pm, dxyz);
                            COO_values.push_back(lapA_val);
                            COO_ROW.push_back(ind);
                            COO_COL.push_back(findex(ix, iy, iz, im, mesh));
                            m_self -= lapA_val;
                        };

                        push_val(i0 + 1, i1, i2, mesh.dx); // +x: ix, ix+1
                        push_val(i0 - 1, i1, i2, mesh.dx); // -x: ix, ix-1
                        push_val(i0, i1 + 1, i2, mesh.dy); // +y: iy, iy+1
                        push_val(i0, i1 - 1, i2, mesh.dy); // -y: iy, iy-1
                        push_val(i0, i1, i2 + 1, mesh.dz); // +z: iz, iz+1
                        push_val(i0, i1, i2 - 1, mesh.dz); // -z: iz, iz-1

                        // add m_i element
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
               Info(), time, time_convert,
               static_cast<double>(af::sparseGetNNZ(matr_CSR)) / static_cast<double>(matr_CSR.elements()));
        fflush(stdout);
    }
    return matr_CSR;
}

af::array PBC_COO_matrix(const double A_exchange, const Mesh& mesh, const bool verbose) {
    return PBC_COO_matrix(std::nullopt, A_exchange, mesh, verbose);
}

af::array PBC_CSR_matrix(const double A_exchange, const Mesh& mesh, const bool verbose) {
    return PBC_CSR_matrix(std::nullopt, A_exchange, mesh, verbose);
}

ExchangeFieldPBC::ExchangeFieldPBC(double A_exchange, Mesh mesh, bool verbose, bool COO)
    : matr(COO ? PBC_COO_matrix(A_exchange, mesh, verbose) : PBC_CSR_matrix(A_exchange, mesh, verbose)) {}

ExchangeFieldPBC::ExchangeFieldPBC(const af::array& A_exchange_field, Mesh mesh, bool verbose, bool COO)
    : matr(COO ? PBC_COO_matrix(A_exchange_field.as(f64), std::nullopt, mesh, verbose)
               : PBC_CSR_matrix(A_exchange_field.as(f64), std::nullopt, mesh, verbose)) {}

// For wrapping only: constructor version taking A_exchange_field
ExchangeFieldPBC::ExchangeFieldPBC(long int A_exchange_field_ptr, Mesh mesh, bool verbose)
    : matr(PBC_CSR_matrix((*(new af::array(*((void**)A_exchange_field_ptr)))).as(f64), std::nullopt, mesh, verbose)) {}

} // namespace magnumafcpp
