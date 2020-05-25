#pragma once
#include "../state.hpp"
#include "LLGTerm.hpp"
#include "arrayfire.h"

namespace magnumafcpp {

class SparseExchangeField : public LLGTerm {
  public:
    SparseExchangeField(double A_exchange, Mesh, bool verbose = true,
                        bool COO = true);
    SparseExchangeField(const af::array& A_exchange_field, Mesh,
                        bool verbose = true, bool COO = true);
    SparseExchangeField(long int A_exchange_field_ptr, Mesh mesh,
                        bool verbose = true);

    af::array h(const State& state); // Field contribution

    double get_cpu_time() { return af_time; } // af time

    const af::array matr;

  private:
    af::array calc_CSR_matrix(const double A_exchange, const Mesh&,
                              const bool verbose);
    af::array calc_COO_matrix(const double A_exchange, const Mesh&,
                              const bool verbose);
    af::array calc_CSR_matrix(const af::array& A_exchange_field, const Mesh&,
                              const bool verbose);
    af::array calc_COO_matrix(const af::array& A_exchange_field, const Mesh&,
                              const bool verbose);
    //int findex(int i0, int i1, int i2, int im, Mesh mesh);
    int findex(uint32_t i0, uint32_t i1, uint32_t i2, uint32_t im,
               const Mesh& mesh);
    double af_time{0};
};
} // namespace magnumafcpp
