#ifndef DEMAG_H
#define DEMAG_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "mesh.hpp"
#include "param.hpp"
#include "func.hpp"

class DemagSolver : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double* get_cpu_time(){return &cpu_time;}

    Param param;
    Mesh mesh;

    DemagSolver (Mesh, Param);
    af::array Nfft;
    af::array h_field;
    af::array hfft;
    af::array mfft;

    double cpu_time{0.};
    af::timer timer_demagsolve;
};
#endif
