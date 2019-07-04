#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

class ExchangeField : public LLGTerm {
  public:
    ExchangeField (double A);
    ExchangeField (af::array A_field);
    ExchangeField (long int A_field_ptr);
    //Field contribution
    af::array h(const State& state);
    //Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field

    double get_cpu_time(){return computation_time_heff;}//!< accumulated heff computation time in [s]

    const double A{0}; //!< Exchange energy in [J/m]
    const af::array A_field{af::array()};

    double computation_time_heff{0.};
    af::timer timer_exchsolve;
    double     time_conv{0.};
    af::timer timer_conv;
    double     time_edges{0.};
    af::timer timer_edges;
};
