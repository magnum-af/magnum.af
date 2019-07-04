#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"


class DmiField : public LLGTerm {
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
    DmiField (Mesh, Material);
    af::array filtr_fd1;
    void correct_edges(af::array& out, const af::array& in);
    af::array n;
    double     cpu_time{0.};
    af::timer timer_dmi;
};
