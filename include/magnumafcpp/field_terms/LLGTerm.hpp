#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "state.hpp"
#include <memory>

namespace magnumafcpp {

// Abstract basis class for all terms in the LLG equation.
class LLGTerm {
  public:
    virtual ~LLGTerm() = default;
    virtual af::array h(const State& state) const = 0;
    /// Calculating the micromagnetic energy \f$E\f$.
    virtual double E(const State& state) const = 0;
    ///< Calculating the micromagnetic energy for a already calculated h field
    ///< (to save computational cost)
    virtual double E(const State& state, const af::array& h) const = 0;
    double get_cpu_time() const { return accumulated_time; };

    /// For wrapping only: pointer to h()
    virtual long int h_ptr(const State& state) { return (long int)(new af::array(h(state)))->get(); }

  protected:
    mutable double accumulated_time{0.};
};

// Aliases used to initialize objects wich inherit from this class
using LlgTerm = std::unique_ptr<LLGTerm>;
using LlgTerms = std::vector<LlgTerm>;

} // namespace magnumafcpp
