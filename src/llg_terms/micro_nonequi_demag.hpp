#ifndef MICRO_NONEQUI_DEMAG_H
#define MICRO_NONEQUI_DEMAG_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../misc.hpp"
#include "../func.hpp"
#include <iostream>
#include <thread>

class NonEquiDemagField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    Material material;
    Mesh mesh;

    NonEquiDemagField (Mesh, Material, bool verbose = false, bool caching = true, unsigned nthreads = 0);
    ///< Array storing the Fourier transfrom of the demag tensor.
    af::array Nfft;
    af::array todel_N;// TODO todel

    double cpu_time{0.};
    af::timer timer_demagsolve;
    
    //For wrapping
    void print_Nfft();
    private:
        const unsigned nthreads;
        af::array N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz);
};
#endif
