#pragma once
#include "arrayfire.h"
#include "field_term.hpp"
#include "nonequispaced_mesh.hpp"

namespace magnumafcpp {

class NonequiTerm : public FieldTerm {
  public:
    NonequiTerm(NonequiMesh nemesh) : nemesh(nemesh){};
    virtual ~NonequiTerm() = default;

    /// Energy calculation: Edemag = - mu0/2 * integral(M . Hdemag) dx
    /// Calculate nonequi distant mesh integral:  integral(M * Hdemag) dx, where
    /// M = Ms * m
    virtual double impl_E_in_J(const State& state, const af::array& h) const override {
        return integral_nonequimesh(h * state.m, state);
    }

  protected:
    NonequiMesh nemesh;

  private:
    double integral_nonequimesh(const af::array& h_times_m,

                                const State& state) const {
        af::array z_spacing_afarray = af::array(1, 1, nemesh.nz, 1, nemesh.z_spacing.data());
        af::array ms_h_times_m;

        // Global or local Ms switch
        if (state.Ms_field.isempty() == true) {
            ms_h_times_m = state.Ms * h_times_m;
            ;
        } else {
            ms_h_times_m = state.Ms_field * h_times_m;
        }

        af::array xy_integral = af::sum(af::sum(af::sum(ms_h_times_m, 0), 1), 3) * nemesh.dx * nemesh.dy;
        af::array xyz_integral = af::sum(xy_integral * z_spacing_afarray, 2);
        return -constants::mu0 / 2. * xyz_integral.as(f64).scalar<double>();
    }
};
} // namespace magnumafcpp
