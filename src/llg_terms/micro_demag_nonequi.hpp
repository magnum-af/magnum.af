#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumaf{


class NonEquiDemagField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    float E(const State& state);
    float E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    float get_cpu_time(){return cpu_time;}

    NonEquiDemagField (NonequispacedMesh nonequimesh, bool verbose = true, bool caching = false, unsigned nthreads = 0);

    af::array Nfft;//!< Array storing the Fourier transfrom of the demag tensor.

    float cpu_time{0.};
    af::timer timer_demagsolve;

    //For wrapping
    void print_Nfft();

    //std::vector<float> z_spacing;
    private:
        const unsigned nthreads;
        af::array calculate_N(int n0_exp, int n1_exp, int n2_exp, float dx, float dy, const std::vector<float> z_spacing);
};
}// namespace magnumaf
