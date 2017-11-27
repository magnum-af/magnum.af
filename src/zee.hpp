#ifndef ZEE_H
#define ZEE_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "state.hpp"
#include "func.hpp"

class Zee : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double get_cpu_time(){return cpu_time;}

    Zee(af::array zee_in, Mesh mesh_in, Param param_in);
    Zee(long int zee_in_addr, Mesh mesh_in, Param param_in);
    af::array zee_field;
    Mesh mesh;
    Param param;
    double cpu_time{0.};
    af::timer timer;
};


#endif
