#pragma once
#include "../state.hpp"
#include "../util/named_type.hpp"
#include "arrayfire.h"
#include "llg_term_mesh_base.hpp"

namespace magnumafcpp {

using RKKY_values = NamedType<const af::array&, struct NamedTypeRKKY_values>;
using Exchange_values =
    NamedType<const af::array&, struct NamedTypeRKKY_values>;

/// Combined Field for RKKY interaction between z-neighbour-cells and exchange
/// interaction between not-RKKY-coupled cells.

///
/// RKKY_values and Exchange_values expect an array with RKKY and exchange
/// constants for each cell, respectively. Both material parameters can be zero
/// or non-zero, where zero values can be used to decouple certain cells. If of
/// two z-neighbouring cells have nonzero RKKY constants they are coupled by \f[
/// H_{RKKY} = \frac{2}{\mu_0 M_s} \sum_i \frac{C_i(\boldsymbol{m}_i -
/// \boldsymbol{m})}{\Delta_i^2} \f], where C_i are the RKKY constants. All
/// other cells are exchange coupled using the values given in the
/// Exchange_values array, provided they have non-zero constants. Jumps in the
/// exchange constants or the RKKY constants are both treated by a harmonic
/// mean.
// TODO this currently does not include jumping M_s as in [mumax3, eqn. (9)],
// but multiplies M_s at the end and hanldes zero values.
///
/// Example: mesh of [n, n, 4]
/// RKKY_values= [0, c_r, c_r, 0] for all n
/// exch_values= [c_ex, c_ex, c_ex, c_ex] for all n
/// layers 1 and 2 are RKKY coupled (counting from 0 to 3)
/// layers 0 and 1 as well as 2 and 3 are exchange coupled
/// All cells are exchange coupled along the xy plane.
/// @param COO Switch between COO and CSR sparse matrix format, true is COO.

class RKKYExchangeField : public LLGTermMeshBase {
  public:
    RKKYExchangeField(RKKY_values rkky_values, Exchange_values exchange_values,
                      Mesh mesh, const af::array& rkky_indices = af::array(),
                      bool verbose = true, bool COO = true);
    RKKYExchangeField(long int rkky_values, long int exchange_values, Mesh mesh,
                      long int rkky_indices, bool verbose = true);

    af::array h(const State& state); // Field contribution

    double get_cpu_time() { return af_time; } // af time

    const af::array matr; // TODO private
  private:
    af::array calc_CSR_matrix(const af::array& RKKY_field,
                              const af::array& A_exchange_field, const Mesh&,
                              const af::array& rkky_indices,
                              const bool verbose);
    af::array calc_COO_matrix(const af::array& RKKY_field,
                              const af::array& A_exchange_field, const Mesh&,
                              const af::array& rkky_indices,
                              const bool verbose);
    int findex(unsigned i0, unsigned i1, unsigned i2, unsigned im,
               const Mesh& mesh);
    double af_time{0};
};
} // namespace magnumafcpp
