#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

namespace magnumaf{



class NonequiUniaxialAnisotropyField : public LLGTerm {
  public:
    NonequiUniaxialAnisotropyField (float Ku1, std::array<float, 3> Ku1_axis = {0, 0, 1});
    NonequiUniaxialAnisotropyField (af::array Ku1_field, std::array<float, 3> Ku1_axis = {0, 0, 1});
    NonequiUniaxialAnisotropyField (af::array Ku1_field, af::array Ku1_axis_field);
    NonequiUniaxialAnisotropyField (float Ku1, float Ku1_axis_0, float Ku1_axis_1, float Ku1_axis_2);
    NonequiUniaxialAnisotropyField (long int Ku1, float Ku1_axis_0, float Ku1_axis_1, float Ku1_axis_2);//!< wrapping only
    NonequiUniaxialAnisotropyField (long int Ku1_field_ptr, long int Ku1_axis_field_ptr);//!< wrapping only
    NonequiUniaxialAnisotropyField (float Ku1, long int Ku1_axis_field_ptr);//!< wrapping only

    af::array h(const State& state);//Field contribution
    long int h_ptr(const State& state);// For wrapping only: pointer to heff
    float E(const State& state);//Energy contribution
    float E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field
    float get_cpu_time(){return computation_time_heff;}//!< accumulated heff computation time in [s]

    const float Ku1{0};//!< [J/m^3]  Uniaxial Anisotropy

    const af::array Ku1_field; //!< Spacially varying anisotropy energy in [J/m^3] defined at each node
    long int get_Ku1_field();

    const std::array<float, 3> Ku1_axis = {0, 0, 0};//!< Anisotropy axis
    const af::array Ku1_axis_field; //!< Spacially varying anisotropy axis
    float get_ku1_axis(int i);//For wrapping only

  private:
    float computation_time_heff{0.};
    std::array<float, 3> get_normalized_vector(std::array<float, 3> vector);
    af::array calc_heff(const State& state);
};
}// namespace magnumaf
