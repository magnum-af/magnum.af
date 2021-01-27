#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "math.hpp"
#include "state.hpp"

namespace magnumafcpp {

/// Combined bulk DMI and exchange Field. Implements
/// \f[ H_{eff} = \frac{-1}{\mu_0 M_s} (2 D_b \ \nabla \times \boldsymbol{m} -2 \nabla \cdot (A \nabla \boldsymbol{m})
/// \f], where \f$D_b\f$ is the bulk DMI constant in [J/m2] and \f$A\f$ is the exchange constant in [J/m].
/// This term is combined as the boundary conditions are coupled:
/// \f[ \boldsymbol{m} \times (2A \frac{\partial
/// \boldsymbol{m}}{\partial \boldsymbol{n}} - D_b \boldsymbol{n} \times \boldsymbol{m}) = 0 \f]
/// , which can be rewritten as
/// \f[ \frac{\partial \boldsymbol{m}}{\partial \boldsymbol{n}} = \frac{D_b}{2A}
/// ( \boldsymbol{n} \times \boldsymbol{m} )\f]
class BulkDMIExchangeField : public MicroTerm {
  public:
    /// \param D_bulk Bulk DMI constant in [J/m2]
    /// \param A Exchange constant in [J/m]
    BulkDMIExchangeField(double D_bulk, double A);

  private:
    double D_bulk_;
    double A_;
    virtual af::array impl_H_in_Apm(const State& state) const override;
};

} // namespace magnumafcpp
