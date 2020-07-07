#pragma once
#include "LLGTerm.hpp"
#include "arrayfire.h"
#include "nonequispaced_mesh.hpp"

namespace magnumafcpp {

class NonequiTermBase : public LLGTerm {
  public:
    virtual ~NonequiTermBase(){};

    // Energy calculation: Edemag = - mu0/2 * integral(M . Hdemag) dx
    // Calculate nonequi distant mesh integral:  integral(M * Hdemag) dx, where
    // M = Ms * m
    double E(const State& state) {
        return integral_nonequimesh(h(state) * state.m, state.nonequimesh,
                                    state);
    }

    double E(const State& state, const af::array& h) {
        return integral_nonequimesh(h * state.m, state.nonequimesh, state);
    }

  private:
    double integral_nonequimesh(const af::array& h_times_m,
                                const NonequispacedMesh& nonequimesh,
                                const State& state) const {
        af::array z_spacing_afarray =
            af::array(1, 1, nonequimesh.nz, 1, nonequimesh.z_spacing.data());
        af::array ms_h_times_m;

        // Global or local Ms switch
        if (state.Ms_field.isempty() == true) {
            ms_h_times_m = state.Ms * h_times_m;
            ;
        } else {
            ms_h_times_m = state.Ms_field * h_times_m;
        }

        af::array xy_integral =
            af::sum(af::sum(af::sum(ms_h_times_m, 0), 1), 3) * nonequimesh.dx *
            nonequimesh.dy;
        af::array xyz_integral = af::sum(xy_integral * z_spacing_afarray, 2);
        return -constants::mu0 / 2. * xyz_integral.scalar<double>();
    }
};
} // namespace magnumafcpp
