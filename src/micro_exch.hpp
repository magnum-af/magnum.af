#ifndef MICRO_EXCH_H
#define MICRO_EXCH_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "state.hpp"
#include "func.hpp"

class ExchSolver : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Field contribution with edges for Energy calculation
    af::array h_withedges(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double* get_cpu_time(){return &cpu_time;}


    Param param;
    Mesh mesh;
    ExchSolver (Mesh, Param);
    af::array filtr;

    double     cpu_time{0.};
    af::timer timer_exchsolve;
    double     time_conv{0.};
    af::timer timer_conv;
    double     time_edges{0.};
    af::timer timer_edges;

    af::array matr;
    int findex(int i0, int i1, int i2, int im, int id);
    int findex(int i0, int i1, int i2, int im);
};
#endif
