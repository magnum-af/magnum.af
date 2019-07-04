#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../nonequispaced_mesh.hpp"
#include "arrayfire.h"

/// Nonequispaced Exchange Field.

/// The second derivative for a nonequispaced mesh can be written as [1, eq. 1.4]: \f[ f_i^{''} = \frac{2 \Big[f_{i-1} + \frac{h_i}{h_{i-1}} f_{i-1} - \big(1 + \frac{h_i}{h_{i-1}} f_i\big)\Big]}{h_i h_{i-1} \big(1 + \frac{h_i}{h_{i-1}}\big)} - \frac{h_i - h_{i-1}}{3} f_i^{'''} + ...\f]
///
/// We calculate the second derivative skipping the term containing f_i as it drops out in the cross product of the LLG: \f[ f_i^{''} = \frac{2 \Big[f_{i-1} + \frac{h_i}{h_{i-1}} f_{i-1}\Big]}{h_i h_{i-1} \big(1 + \frac{h_i}{h_{i-1}}\big)}\f]
///
/// [1] "A simple finite-difference grid with non-constant intervals" by Hilding Sundqvist & George Veronis
class NonequiExchangeField : public LLGTerm {
  public:
    NonequiExchangeField (double A_exchange, NonequispacedMesh, bool verbose = true);
    NonequiExchangeField (const af::array& A_exchange_field, NonequispacedMesh, bool verbose = true);
    NonequiExchangeField (long int A_exchange_field_ptr, NonequispacedMesh mesh, bool verbose = true);

    af::array h(const State& state);//Field contribution
    double E(const State& state);//Energy contribution
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field

    double get_cpu_time(){return af_time;}//af time

    //const double A_exchange { 0 };
    //const af::array A_exchange_field = af::array();// empty array if not specified in constructor


  private:
    const af::array matr;
    af::array calc_CSR_matrix(const double A_exchange, const NonequispacedMesh&, const bool verbose);
    af::array calc_CSR_matrix(const af::array& A_exchange_field, const NonequispacedMesh&, const bool verbose);
    int findex(int i0, int i1, int i2, int im, NonequispacedMesh mesh);
    double af_time { 0 };
};
