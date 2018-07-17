#ifndef ATOMISTIC_DMI_H
#define ATOMISTIC_DMI_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
class ATOMISTIC_DMI : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    //CPU time
    double get_cpu_time(){return cpu_time;}

    //ATOMISTIC_DMI ();
    ATOMISTIC_DMI (const Mesh& mesh, const Param& param);
    af::array n;
    af::array filtr_fd1;
    //af::array filtr_atom_dmi_x;
    //af::array filtr_atom_dmi_y;
    //af::array filtr_atom_dmi_z;

    double     cpu_time{0.};
    af::timer timer_dmi;
};
#endif
