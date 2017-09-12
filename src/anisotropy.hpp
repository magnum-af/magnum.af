#ifndef ANISOTROPY_H
#define ANISOTROPY_H
#include "arrayfire.h"
#include "mesh.hpp"
#include "param.hpp"
#include "func.hpp"
#include "LLGTerm.hpp"
class ANISOTROPY : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double* get_cpu_time(){return &cpu_time;}

    ANISOTROPY (Mesh, Param);
    Param param;
    Mesh mesh;

    af::array eu;//Uniaxial anisotropy normal vector

    double     cpu_time{0.};
    af::timer timer_anisotropy;
};
#endif
