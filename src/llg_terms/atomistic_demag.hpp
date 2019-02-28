#ifndef AtomisticDipoleDipoleField_H
#define AtomisticDipoleDipoleField_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"

class AtomisticDipoleDipoleField: public LLGTerm{
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    AtomisticDipoleDipoleField (Mesh);

    af::array Nfft;
    double cpu_time{0.};
    af::timer timer_demagsolve;
};
#endif
