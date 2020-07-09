#pragma once
#include "LLGTerm.hpp"
#include "arrayfire.h"
#include "nonequispaced_mesh.hpp"

namespace magnumafcpp {

class NonequiTermBase : public LLGTerm {
  public:
    NonequiTermBase(NonequispacedMesh nemesh) : nemesh(nemesh){};
    virtual ~NonequiTermBase(){};

    // Energy calculation: Edemag = - mu0/2 * integral(M . Hdemag) dx
    // Calculate nonequi distant mesh integral:  integral(M * Hdemag) dx, where
    // M = Ms * m
    double E(const State& state) {
        return integral_nonequimesh(h(state) * state.m, state);
    }

    double E(const State& state, const af::array& h) {
        return integral_nonequimesh(h * state.m, state);
    }

  private:
    const NonequispacedMesh nemesh;
    double integral_nonequimesh(const af::array& h_times_m,

                                const State& state) const {
        af::array z_spacing_afarray =
            af::array(1, 1, nemesh.nz, 1, nemesh.z_spacing.data());
        af::array ms_h_times_m;

        // Global or local Ms switch
        if (state.Ms_field.isempty() == true) {
            ms_h_times_m = state.Ms * h_times_m;
            ;
        } else {
            ms_h_times_m = state.Ms_field * h_times_m;
        }

        af::array xy_integral =
            af::sum(af::sum(af::sum(ms_h_times_m, 0), 1), 3) * nemesh.dx *
            nemesh.dy;
        af::array xyz_integral = af::sum(xy_integral * z_spacing_afarray, 2);
        return -constants::mu0 / 2. * xyz_integral.scalar<double>();
    }
};
} // namespace magnumafcpp
