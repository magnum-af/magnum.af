#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumaf{



class DemagField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //CPU time
    float get_cpu_time(){return cpu_time;}

    DemagField (Mesh, bool verbose = false, bool caching = true, unsigned nthreads = 0);
    ///< Array storing the Fourier transfrom of the demag tensor.
    af::array Nfft;

    float cpu_time{0.};
    af::timer timer_demagsolve;

    //For wrapping
    void print_Nfft();
    private:
        const unsigned nthreads;
        af::array N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, float dx, float dy, float dz);
};
}// namespace magnumaf
