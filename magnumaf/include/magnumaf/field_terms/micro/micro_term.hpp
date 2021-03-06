#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "field_term.hpp"
#include "state.hpp"
#include <memory>

namespace magnumaf {

class MicroTerm : public FieldTerm {
  public:
    virtual ~MicroTerm() = default;

  private:
    /// Calculating the micromagnetic energy \f$E\f$.
    /// This is a prototype for all llgterms with are linear in m and must be
    /// overwritten in e.g. zeeman where factor 1/2 becomes 1.
    virtual double impl_E_in_J(const State& state, const af::array& h) const override {
        if (state.Ms_field.isempty()) {
            return -constants::mu0 / 2. * state.Ms *
                   af::sum(af::sum(af::sum(af::sum(h * state.m, 0), 1), 2), 3).as(f64).scalar<double>() *
                   state.mesh.dx * state.mesh.dy * state.mesh.dz;
        } else {
            return -constants::mu0 / 2. *
                   af::sum(af::sum(af::sum(af::sum(state.get_Ms_field_in_vec_dims() * h * state.m, 0), 1), 2), 3)
                       .as(f64)
                       .scalar<double>() *
                   state.mesh.dx * state.mesh.dy * state.mesh.dz;
        }
    }
};

} // namespace magnumaf
