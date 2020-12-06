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
    ///< Calculating the micromagnetic energy from the h field
    virtual double E(const State& state, const af::array& h) const = 0;
    /// Calculating the micromagnetic energy \f$E\f$.
    // virtual double E(const State& state) const = 0;
    double E(const State& state) const { return E(state, h(state)); };
    auto h_and_E(const State& state) {
        const auto htmp = h(state);
        return std::make_pair(htmp, E(state, htmp));
    };

    double get_cpu_time() const { return accumulated_time; };

    /// For wrapping only: pointer to h()
    virtual long int h_ptr(const State& state) { return (long int)(new af::array(h(state)))->get(); }

  protected:
    mutable double accumulated_time{0.};
};

// helper function to create unique_ptrs
template <typename T> std::unique_ptr<LLGTerm> cp_to_uptr(const T& t) { return std::make_unique<T>(t); }

/// returns std::vector<std::unique_ptr<LLGTerm>> from args
/// args called with std::move() are moved, copied elsewise
template <typename... Args> auto to_vec(Args... args) {
    std::vector<std::unique_ptr<LLGTerm>> v;
    (v.push_back(std::unique_ptr<LLGTerm>(std::make_unique<Args>(args))), ...);
    return v;
}

// Always moves arguments (except when arg is const, then copies)
template <typename T> std::unique_ptr<LLGTerm> to_uptr(T t) { return std::make_unique<T>(t); }
template <typename... Args> auto mv_to_vec(Args&&... args) {
    std::vector<std::unique_ptr<LLGTerm>> v;
    (v.push_back(to_uptr(std::move(args))), ...); // works and moves as well
    return v;
}

// Aliases used to initialize objects wich inherit from this class
using LlgTerm = std::unique_ptr<LLGTerm>;
using LlgTerms = std::vector<LlgTerm>;

} // namespace magnumafcpp
