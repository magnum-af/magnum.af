#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumaf {

class SparseExchangeField : public MicroTerm {
  public:
    SparseExchangeField(double A_exchange, Mesh, bool verbose = true, bool COO = true);
    SparseExchangeField(const af::array& A_exchange_field, Mesh, bool verbose = true, bool COO = true);
    SparseExchangeField(long int A_exchange_field_ptr, Mesh mesh, bool verbose = true);

    af::array get_matr() const { return matr; };
    long int pywrap_get_sparse_matrix_ptr() const { return util::pywrap::send_copy_to_py(this->matr); }

  private:
    af::array matr;
    virtual af::array impl_H_in_Apm(const State& state) const override; // Field contribution
};
} // namespace magnumaf
