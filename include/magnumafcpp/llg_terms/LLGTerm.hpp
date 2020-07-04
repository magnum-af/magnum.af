#pragma once
#include "../constants.hpp"
#include "../state.hpp"
#include "arrayfire.h"
#include <memory>

namespace magnumafcpp {

// Abstract basis class for all terms in the LLG equation.
class LLGTerm {
  public:
    virtual af::array h(const State& state) = 0;
    /// Calculating the micromagnetic energy \f$E\f$.
    virtual double E(const State& state) = 0;
    ///< Calculating the micromagnetic energy for a already calculated h field
    ///< (to save computational cost)
    virtual double E(const State& state, const af::array& h) = 0;
    virtual double get_cpu_time() = 0;

    /// For wrapping only: pointer to h()
    virtual long int h_ptr(const State& state) {
        return (long int)(new af::array(h(state)))->get();
    }
    virtual ~LLGTerm(){};
};

// Aliases used to initialize objects wich inherit from this class
using LlgTerm = std::shared_ptr<LLGTerm>;
using LlgTerms = std::vector<LlgTerm>;

} // namespace magnumafcpp
