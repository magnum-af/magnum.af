#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "state.hpp"
#include <memory>

namespace magnumafcpp {

// Abstract basis class for all terms in the LLG equation.
class Fieldterm {
  public:
    virtual ~Fieldterm() = default;
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

// helper functions:

// create a unique_ptr (e.g. from ctor or copy-ctor)
template <typename T, class... Args> auto fieldterm_ptr(Args&&... args) {
    return std::unique_ptr<Fieldterm>(std::make_unique<T>(std::forward<Args>(args)...));
}

// create unique_ptr to a copy
template <typename T> std::unique_ptr<Fieldterm> cp_to_uptr(const T& t) { return std::make_unique<T>(t); }

// moves to a unique_ptr, do not use element afterwards
// requires std::move() at callsite
template <typename T> std::unique_ptr<Fieldterm> mv_to_uptr(T&& t) { return std::make_unique<T>(t); }

/// returns std::vector<std::unique_ptr<Fieldterm>> from args
/// args called with std::move() are moved, copied elsewise
template <typename... Args> auto to_vec(Args... args) {
    std::vector<std::unique_ptr<Fieldterm>> v;
    (v.push_back(std::unique_ptr<Fieldterm>(std::make_unique<Args>(args))), ...);
    return v;
}

// Always moves arguments, except when arg is const itself
// Uses mv-ctor on args or copy-ctor as fallback (when arg is const)
template <typename T> std::unique_ptr<Fieldterm> to_uptr(T t) { return std::make_unique<T>(t); }
template <typename... Args> auto mv_to_vec(Args&&... args) {
    std::vector<std::unique_ptr<Fieldterm>> v;
    (v.push_back(to_uptr(std::move(args))), ...);
    return v;
}

// Aliases used to initialize objects wich inherit from this class
using uptr_Fieldterm = std::unique_ptr<Fieldterm>;
using vec_uptr_Fieldterm = std::vector<uptr_Fieldterm>;

} // namespace magnumafcpp
