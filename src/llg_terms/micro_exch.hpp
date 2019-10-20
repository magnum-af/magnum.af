#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumaf{


class ExchangeField : public LLGTerm {
  public:
    ExchangeField (float A);
    ExchangeField (af::array A_field);
    ExchangeField (long int A_field_ptr);
    //Field contribution
    af::array h(const State& state);
    //Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state);
    //Energy contribution
    float E(const State& state);

    float get_cpu_time(){return computation_time_heff;}//!< accumulated heff computation time in [s]

    const float A{0}; //!< Exchange energy in [J/m]
    const af::array A_field{af::array()};

    float computation_time_heff{0.};
    af::timer timer_exchsolve;
    float     time_conv{0.};
    af::timer timer_conv;
    float     time_edges{0.};
    af::timer timer_edges;
};
}// namespace magnumaf
