#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumaf {

///
/// Slonczewski Spin Transfer Torque model
/// \f[ \boldsymbol{H}_{stt} = - \frac{j_e \hbar}{2 e \mu_0 M_s t_{fl}}
/// [\eta_{damp} \boldsymbol{m} \times \boldsymbol{p} + \eta_{field} \boldsymbol{p}] \f]
/// where \f$j_e\f$ is the current in [A], \f$\hbar\f$ Planck constant, \f$e\f$ is the electron charge,
/// \f$\eta_{damp}\f$ and \f$\eta_{field}\f$ are damping and field like constants in [a.u.], \f$\boldsymbol{p}\f$ is the
/// spin polarization as unit length vector field in [a.u.], \f$t_{fl}\f$ is the free layer thickness in [m].
///

class SpinTransferTorqueField : public MicroTerm {
  public:
    /// \param polarization unit-less polarization vector field in []
    /// \param eta_dampinglike unit-less damping like factor in []
    /// \param eta_field unit-less field like factor in []
    /// \param j_e current in [A]
    /// \param fl_thickness free-layer thickness in [m]
    SpinTransferTorqueField(const af::array& polarization, double eta_damping, double eta_field, double j_e,
                            double fl_thickness);
    SpinTransferTorqueField(long int polarization_ptr, double eta_damping, double eta_field, double j_e,
                            double fl_thickness);

    af::array polarization; ///< spin polarization as cell-dependent unit length vector field in [a.u.]
    double eta_damping;     ///< Damping like constant in [a.u.]
    double eta_field;       ///< Field like constant in [a.u.]
    double j_e;             ///< Current in [A]
    double fl_thickness;    ///< Free layer thickness in [m]

  private:
    virtual af::array impl_H_in_Apm(const State& state) const override; ///< Effective field term contribution
};
} // namespace magnumaf
