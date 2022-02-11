#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumaf {

///
/// Zhang Li Spin Transfer Torque model
///

class SpinTransferTorqueZhangLiField : public MicroTerm {
  public:
    /// \param j : cell dependent spin polarization vector of dims [nx ny nz 3] in [A/m^2]
    /// \param beta : dimension-less polarization rate
    /// \param xi : degree of nonadiabacity
    SpinTransferTorqueZhangLiField(const af::array& j, double beta, double xi) : j_(j), beta_(beta), xi_(xi) {}
    SpinTransferTorqueZhangLiField(long int j_ptr, double beta, double xi);

    af::array j_; ///< spin polarization as cell-dependent vector field in [A/m^2]
    double beta_; ///< Polarization rate of the conducting electrons in [a.u.]
    double xi_;   ///< Degree of nonadiabacity in [a.u.]

  private:
    virtual af::array impl_H_in_Apm(const State& state) const override;
};
} // namespace magnumaf
