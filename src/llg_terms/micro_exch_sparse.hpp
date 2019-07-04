#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"

class SparseExchangeField : public LLGTerm {
  public:
    SparseExchangeField (double A_exchange, Mesh, bool verbose = true);
    SparseExchangeField (const af::array& A_exchange_field, Mesh, bool verbose = true);
    SparseExchangeField (long int A_exchange_field_ptr, Mesh mesh, bool verbose = true);

    af::array h(const State& state);//Field contribution
    double E(const State& state);//Energy contribution
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field

    double get_cpu_time(){return af_time;}//af time

    //const double A_exchange { 0 };
    //const af::array A_exchange_field = af::array();// empty array if not specified in constructor


  private:
    const af::array matr;
    af::array calc_CSR_matrix(const double A_exchange, const Mesh&, const bool verbose);
    af::array calc_CSR_matrix(const af::array& A_exchange_field, const Mesh&, const bool verbose);
    int findex(int i0, int i1, int i2, int im, Mesh mesh);
    double af_time { 0 };
};
