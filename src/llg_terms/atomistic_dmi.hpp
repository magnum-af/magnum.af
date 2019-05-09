#ifndef AtomisticDmiField_H
#define AtomisticDmiField_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
class AtomisticDmiField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    //AtomisticDmiField ();
    AtomisticDmiField (const Mesh& mesh, const Material& material);
    af::array n;
    af::array filtr_fd1;
    //af::array filtr_atom_dmi_x;
    //af::array filtr_atom_dmi_y;
    //af::array filtr_atom_dmi_z;

    double     cpu_time{0.};
    af::timer timer_dmi;
};
#endif
