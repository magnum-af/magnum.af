#ifndef SPARSE_EXCHANGE_FIELD_H
#define SPARSE_EXCHANGE_FIELD_H
#include "arrayfire.h"
#include "LLGTerm.hpp"
#include "../state.hpp"
#include "../func.hpp"

class SparseExchangeField : public LLGTerm {
  public:
    SparseExchangeField (Mesh, Material);

    af::array h(const State& state);//Field contribution
    double E(const State& state);//Energy contribution
    double E(const State& state, const af::array& h);///< Calculating the micromagnetic energy for a already calculated h field

    double get_cpu_time(){return af_time;}//af time


  private:
    af::array matr;
    int findex(int i0, int i1, int i2, int im, Mesh mesh);
    double af_time{0.};
};
#endif
