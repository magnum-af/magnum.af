#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
namespace magnumaf {
class DemagFieldPBC : public MicroTerm {
  private:
    virtual af::array impl_H_in_Apm(const State& state) const override;
};
} // namespace magnumaf
