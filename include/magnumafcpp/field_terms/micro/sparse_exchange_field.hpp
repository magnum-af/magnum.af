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

    virtual af::array impl_H_in_Apm(const State& state) const override; // Field contribution

    af::array get_matr() const { return matr; };

  private:
    af::array matr;
};
} // namespace magnumafcpp
