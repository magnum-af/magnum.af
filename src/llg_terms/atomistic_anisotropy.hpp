#ifndef AtomisticUniaxialAnisotropyField_H
#define AtomisticUniaxialAnisotropyField_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
class AtomisticUniaxialAnisotropyField : public LLGTerm {
  public:
    //Field contribution
    af::array h(const State& state);
    //Energy contribution
    double E(const State& state);
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    //CPU time
    double get_cpu_time(){return cpu_time;}

    AtomisticUniaxialAnisotropyField (const Mesh& mesh, const Material& material);
    //AtomisticUniaxialAnisotropyField (const Mesh&, const Material&);

    af::array eu;//Uniaxial anisotropy normal vector

    double     cpu_time{0.};
    af::timer timer_anisotropy;
};
#endif
