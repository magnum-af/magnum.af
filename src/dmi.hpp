#ifndef DMI_H
#define DMI_H
#include "arrayfire.h"
#include "mesh.hpp"
#include "param.hpp"
#include "func.hpp"
#include "LLGTerm.hpp"
class DMI : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double* get_cpu_time(){return &cpu_time;}

    Param param;
    Mesh mesh;
    DMI (Mesh, Param);
    af::array filtr_fd1;
    void correct_edges(af::array& out, const af::array& in);
    af::array n;
    double     cpu_time{0.};
    af::timer timer_dmi;
};
#endif
