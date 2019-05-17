#ifndef MICRO_UniaxialAnisotropyField_H
#define MICRO_UniaxialAnisotropyField_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"
class UniaxialAnisotropyField : public LLGTerm {
  public:
    UniaxialAnisotropyField (double Ku1, std::array<double, 3> Ku1_axis = {0,0,1});
    UniaxialAnisotropyField (double Ku1, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2);

    af::array h(const State& state);//Field contribution
    double E(const State& state);//Energy contribution
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    double get_cpu_time(){return computation_time_heff;}//!< accumulated heff computation time in [s]
    const double Ku1;//!< [J/m^3]  Uniaxial Anisotropy 
    const std::array<double, 3> Ku1_axis;//!< Anisotropy axis
    double get_ku1_axis(int i);//For wrapping only

  private:
    af::array eu;//Uniaxial anisotropy normal vector
    double computation_time_heff{0.};
    std::array<double, 3> get_normalized_vector(std::array<double, 3> vector);
};
#endif
