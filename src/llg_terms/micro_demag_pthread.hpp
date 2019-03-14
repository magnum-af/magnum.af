#ifndef MICRO_DEMAG_PTHREAD_H
#define MICRO_DEMAG_PTHREAD_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
#include <pthread.h>
//#include <thread>

class DemagFieldMultithread : public LLGTerm {
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

    DemagFieldMultithread (Mesh, Material, bool verbose = false, bool caching = true, int nthreads = 8);
    ///< Array storing the Fourier transfrom of the demag tensor.
    af::array Nfft;

    double cpu_time{0.};
    af::timer timer_demagsolve;
    
    //For wrapping
    void print_Nfft();
    private:
        //void*  setup_N(void* arg);
        //void setup_N(Mesh mesh);
        //double newellg(double x, double y, double z);
        //double newellf(double x, double y, double z);
        //double Nxxg(int ix, int iy, int iz, double dx, double dy, double dz);
        //double Nxxf(int ix, int iy, int iz, double dx, double dy, double dz);
        const int nthreads;
        af::array N_cpp_alloc(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz);
};
#endif
