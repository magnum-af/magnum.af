#ifndef ATOMISTIC_ANISOTROPY_H
#define ATOMISTIC_ANISOTROPY_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "state.hpp"
#include "func.hpp"
class ATOMISTIC_ANISOTROPY : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double get_cpu_time(){return cpu_time;}

    ATOMISTIC_ANISOTROPY (const Mesh& mesh, const Param& param);
    //ATOMISTIC_ANISOTROPY (const Mesh&, const Param&);

    af::array eu;//Uniaxial anisotropy normal vector

    double     cpu_time{0.};
    af::timer timer_anisotropy;
};
#endif
