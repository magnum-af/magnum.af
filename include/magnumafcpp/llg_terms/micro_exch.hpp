#pragma once
#include "../state.hpp"
#include "arrayfire.h"
#include "integrator_term_mesh_base.hpp"

namespace magnumafcpp {

class ExchangeField : public IntegratorTermMeshBase {
  public:
    using LLGTerm::E;
    ExchangeField(double A);
    ExchangeField(af::array A_field);
    ExchangeField(long int A_field_ptr);
    // Field contribution
    af::array h(const State& state) override;
    // Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state);
    // Energy contribution
    double E(const State& state) override;

    double get_cpu_time() override { return computation_time_heff; } //!< accumulated heff computation time in [s]

    const double A{0}; //!< Exchange energy in [J/m]
    const af::array A_field;

    double computation_time_heff{0.};
    af::timer timer_exchsolve;
    double time_conv{0.};
    af::timer timer_conv;
    double time_edges{0.};
    af::timer timer_edges;
};
} // namespace magnumafcpp
