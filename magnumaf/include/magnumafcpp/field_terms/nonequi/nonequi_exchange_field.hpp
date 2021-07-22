#pragma once
#include "arrayfire.h"
#include "field_terms/nonequi/nonequi_term.hpp"
#include "nonequispaced_mesh.hpp"
#include "state.hpp"

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
class NonequiExchangeField : public NonequiTerm {
  public:
    NonequiExchangeField(const NonequiMesh& nemesh, double A_exchange, bool verbose = true, bool COO = true);
    NonequiExchangeField(const NonequiMesh& nemesh, const af::array& A_exchange_field, bool verbose = true, bool COO = true);
    NonequiExchangeField(const NonequiMesh& nemesh, long int A_exchange_field_ptr, bool verbose = true, bool COO = true);

    af::array get_matr() const { return matr; };

  private:
    af::array matr;
    virtual af::array impl_H_in_Apm(const State& state) const override; // Field contribution
};

} // namespace magnumafcpp
