#pragma once
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "arrayfire.h"
#include "../util/named_type.hpp"

namespace magnumaf{

using RKKY_values = NamedType<const af::array&, struct NamedTypeRKKY_values>;
using Exchange_values = NamedType<const af::array&, struct NamedTypeRKKY_values>;

class RKKYExchangeField : public LLGTerm {
  public:
    //RKKYExchangeField (double A_exchange, Mesh, bool verbose = true);
    //TODO RKKY matrix takes values for interface 
    //e.g. rkky between layer 1 and 2
    //layer 0 : 0
    //layer 1 : rkkyval
    //layer 2 : rkkyval
    //layer 3 : 0
    RKKYExchangeField (RKKY_values rkky_values, Exchange_values exchange_values, Mesh, bool verbose = true);

    af::array h(const State& state);//Field contribution

    double get_cpu_time(){return af_time;}//af time


  private:
    const af::array matr;
    af::array calc_CSR_matrix(const af::array& RKKY_field, const af::array& A_exchange_field, const Mesh&, const bool verbose);
    int findex(int i0, int i1, int i2, int im, Mesh mesh);
    double af_time { 0 };
};
}// namespace magnumaf
