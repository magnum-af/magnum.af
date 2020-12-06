#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

class UniaxialAnisotropyField : public MicroTerm {
  public:
    UniaxialAnisotropyField(double Ku1, std::array<double, 3> Ku1_axis = {0, 0, 1});
    UniaxialAnisotropyField(af::array Ku1_field, std::array<double, 3> Ku1_axis = {0, 0, 1});
    UniaxialAnisotropyField(af::array Ku1_field, af::array Ku1_axis_field);
    UniaxialAnisotropyField(double Ku1, double Ku1_axis_0, double Ku1_axis_1, double Ku1_axis_2);
    UniaxialAnisotropyField(long int Ku1, double Ku1_axis_0, double Ku1_axis_1,
                            double Ku1_axis_2); //!< wrapping only
    UniaxialAnisotropyField(long int Ku1_field_ptr,
                            long int Ku1_axis_field_ptr); //!< wrapping only
    UniaxialAnisotropyField(double Ku1,
                            long int Ku1_axis_field_ptr); //!< wrapping only

    virtual af::array h(const State& state) const override; // Field contribution

    double Ku1{0}; //!< [J/m^3]  Uniaxial Anisotropy

    af::array Ku1_field; //!< Spacially varying anisotropy energy in [J/m^3] defined
                         //!< at each node. Expects size of [nx, ny, nz, 1], stores as
                         //!< [nx, ny, nz, 3];
    long int get_Ku1_field();

    std::array<double, 3> Ku1_axis = {0, 0, 0}; //!< Anisotropy axis
    af::array Ku1_axis_field;                   //!< Spacially varying anisotropy axis
    double get_ku1_axis(int i);                 // For wrapping only

  private:
};
} // namespace magnumafcpp
