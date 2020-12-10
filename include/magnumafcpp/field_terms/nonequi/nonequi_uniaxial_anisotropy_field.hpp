#pragma once
#include "arrayfire.h"
#include "field_terms/nonequi/nonequi_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

class NonequiUniaxialAnisotropyField : public NonequiTerm {
  public:
    NonequiUniaxialAnisotropyField(NonequiMesh nemesh, double Ku1, std::array<double, 3> Ku1_axis = {0, 0, 1});
    NonequiUniaxialAnisotropyField(NonequiMesh nemesh, af::array Ku1_field, std::array<double, 3> Ku1_axis = {0, 0, 1});
    NonequiUniaxialAnisotropyField(NonequiMesh nemesh, af::array Ku1_field, af::array Ku1_axis_field);
    NonequiUniaxialAnisotropyField(NonequiMesh nemesh, double Ku1, double Ku1_axis_0, double Ku1_axis_1,
                                   double Ku1_axis_2);
    NonequiUniaxialAnisotropyField(NonequiMesh nemesh, long int Ku1, double Ku1_axis_0, double Ku1_axis_1,
                                   double Ku1_axis_2); //!< wrapping only
    NonequiUniaxialAnisotropyField(NonequiMesh nemesh, long int Ku1_field_ptr,
                                   long int Ku1_axis_field_ptr); //!< wrapping only
    NonequiUniaxialAnisotropyField(NonequiMesh nemesh, double Ku1,
                                   long int Ku1_axis_field_ptr); //!< wrapping only

    long int get_Ku1_field() const;

    double Ku1{0};              //!< [J/m^3]  Uniaxial Anisotropy
    double get_ku1_axis(int i) const; // For wrapping only
  private:
    af::array Ku1_field{};                      //!< Spacially varying anisotropy energy in
    std::array<double, 3> Ku1_axis = {0, 0, 0}; //!< Anisotropy axis
    af::array Ku1_axis_field{};                 //!< Spacially varying anisotropy axis
    std::array<double, 3> get_normalized_vector(std::array<double, 3> vector);
    af::array calc_heff(const State& state) const;

    virtual af::array impl_H_in_Apm(const State& state) const override; // Field contribution
                                                                        //!< [J/m^3] defined at each node
};
} // namespace magnumafcpp
