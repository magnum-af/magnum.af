#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumaf{



class DmiField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    float E(const State& state);
    float E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    float get_cpu_time(){return cpu_time;}

    Material material;
    Mesh mesh;
    DmiField (Mesh, Material);
    af::array filtr_fd1;
    void correct_edges(af::array& out, const af::array& in);
    af::array n;
    float     cpu_time{0.};
    af::timer timer_dmi;
};
}// namespace magnumaf
