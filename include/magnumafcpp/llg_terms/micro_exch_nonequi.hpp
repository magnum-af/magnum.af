#pragma once
#include "../nonequispaced_mesh.hpp"
#include "../state.hpp"
#include "LLGTerm.hpp"
#include "arrayfire.h"

namespace magnumafcpp {

/// Nonequispaced Exchange Field.

/// The second derivative for a nonequispaced mesh can be written as [1,
/// eq. 1.4]: \f[ f_i^{''} = \frac{2 \Big[f_{i-1} + \frac{h_i}{h_{i-1}} f_{i-1}
/// - \big(1 + \frac{h_i}{h_{i-1}} f_i\big)\Big]}{h_i h_{i-1} \big(1 +
/// \frac{h_i}{h_{i-1}}\big)} - \frac{h_i - h_{i-1}}{3} f_i^{'''} + ...\f]
///
/// We calculate the second derivative skipping the term containing f_i as it
/// drops out in the cross product of the LLG: \f[ f_i^{''} = \frac{2
/// \Big[f_{i-1} + \frac{h_i}{h_{i-1}} f_{i-1}\Big]}{h_i h_{i-1} \big(1 +
/// \frac{h_i}{h_{i-1}}\big)}\f]
///
/// [1] "A simple finite-difference grid with non-constant intervals" by Hilding
/// Sundqvist & George Veronis
class NonequiExchangeField : public LLGTerm {
  public:
    NonequiExchangeField(double A_exchange, NonequispacedMesh,
                         bool verbose = true, bool COO = true);
    NonequiExchangeField(const af::array& A_exchange_field, NonequispacedMesh,
                         bool verbose = true, bool COO = true);
    NonequiExchangeField(long int A_exchange_field_ptr, NonequispacedMesh mesh,
                         bool verbose = true, bool COO = true);

    af::array h(const State& state); // Field contribution
    double E(const State& state);    // Energy contribution
    double E(const State& state,
             const af::array& h); ///< Calculating the micromagnetic energy for
                                  ///< a already calculated h field

    double get_cpu_time() { return af_time; } // af time

    const af::array matr;

  private:
    af::array calc_CSR_matrix(const double A_exchange, const NonequispacedMesh&,
                              const bool verbose);
    af::array calc_CSR_matrix(const af::array& A_exchange_field,
                              const NonequispacedMesh&, const bool verbose);
    af::array calc_COO_matrix(const double A_exchange, const NonequispacedMesh&,
                              const bool verbose);
    af::array calc_COO_matrix(const af::array& A_exchange_field,
                              const NonequispacedMesh&, const bool verbose);
    unsigned findex(unsigned i0, unsigned i1, unsigned i2, unsigned im,
                    const NonequispacedMesh& mesh);
    double af_time{0};
};

} // namespace magnumafcpp
