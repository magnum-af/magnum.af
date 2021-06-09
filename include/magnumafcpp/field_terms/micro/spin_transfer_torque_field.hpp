#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

///
/// Slonczewski Spin Transfer Torque model
/// \f[ \boldsymbol{H}_{stt} = - \frac{j_e \hbar}{2 e \mu_0 M_s t_{fl}}
/// [\nu_{damp} \boldsymbol{m} \times \boldsymbol{p} + \nu_{field} \boldsymbol{p}] \f]
/// where \f$j_e\f$ is the current in [A], \f$\hbar\f$ Planck constant, \f$e\f$ is the electron charge, \f$\nu_{damp}\f$
/// and \f$\nu_{field}\f$ are damping and field like constants in [a.u.], \f$\boldsymbol{p}\f$ is the spin polarization
/// as unit length vector field in [a.u.], \f$t_{fl}\f$ is the free layer thickness in [m].
///

class SpinTransferTorqueField : public MicroTerm {
  public:
    /// \param polarization unit-less polarization vector field in []
    /// \param nu_dampinglike unit-less damping like factor in []
    /// \param nu_field unit-less field like factor in []
    /// \param j_e current in [A]
    /// \param fl_thickness free-layer thickness in [m]
    SpinTransferTorqueField(af::array polarization, double nu_damping, double nu_field, double j_e,
                            double fl_thickness);
    SpinTransferTorqueField(long int polarization_ptr, double nu_damping, double nu_field, double j_e,
                            double fl_thickness);

    af::array polarization; ///< spin polarization as cell-dependent unit length vector field in [a.u.]
    double nu_damping;      ///< Damping like constant in [a.u.]
    double nu_field;        ///< Field like constant in [a.u.]
    double j_e;             ///< Current in [A]
    double fl_thickness{};  ///< Free layer thickness in [m]

  private:
    virtual af::array impl_H_in_Apm(const State& state) const override; ///< Effective field term contribution
};
} // namespace magnumafcpp
