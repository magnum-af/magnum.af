#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

class SparseExchangeField : public MicroTerm {
  public:
    SparseExchangeField(double A_exchange, Mesh, bool verbose = true, bool COO = true);
    SparseExchangeField(const af::array& A_exchange_field, Mesh, bool verbose = true, bool COO = true);
    SparseExchangeField(long int A_exchange_field_ptr, Mesh mesh, bool verbose = true);

    virtual af::array h(const State& state) const override; // Field contribution

    const af::array matr;

  private:
    af::array calc_CSR_matrix(const double A_exchange, const Mesh&, const bool verbose);
    af::array calc_COO_matrix(const double A_exchange, const Mesh&, const bool verbose);
    af::array calc_CSR_matrix(const af::array& A_exchange_field, const Mesh&, const bool verbose);
    af::array calc_COO_matrix(const af::array& A_exchange_field, const Mesh&, const bool verbose);
    unsigned findex(unsigned i0, unsigned i1, unsigned i2, unsigned im, const Mesh& mesh);
    double af_time{0};
};
} // namespace magnumafcpp
