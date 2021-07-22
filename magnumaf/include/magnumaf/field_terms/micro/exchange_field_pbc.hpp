#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumaf {

class ExchangeFieldPBC : public MicroTerm {
  public:
    ExchangeFieldPBC(double A_exchange, Mesh, bool verbose = true, bool COO = true);
    ExchangeFieldPBC(const af::array& A_exchange_field, Mesh, bool verbose = true, bool COO = true);
    ExchangeFieldPBC(long int A_exchange_field_ptr, Mesh mesh, bool verbose = true);
    af::array get_matr() const { return matr; };

  private:
    af::array matr;
    virtual af::array impl_H_in_Apm(const State& state) const override; // Field contribution
};
} // namespace magnumaf
