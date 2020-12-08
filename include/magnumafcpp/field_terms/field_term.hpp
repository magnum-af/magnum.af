#pragma once
#include "arrayfire.h"
#include "constants.hpp"
#include "state.hpp"
#include <memory>
#include <numeric>

namespace magnumafcpp {

// Abstract basis class for all terms in the LLG equation.
class FieldTerm {
  public:
    virtual ~FieldTerm() = default;
    af::array H_eff(const State& state) const;
    double Energy(const State& state, const af::array& h) const;
    double Energy(const State& state) const;
    /// Calculating the micromagnetic energy \f$E\f$.
    // virtual double E(const State& state) const = 0;

    std::pair<af::array, double> h_and_E(const State& state) const {
        const auto htmp = H_eff(state);
        return {htmp, Energy(state, htmp)};
    };

    double get_cpu_time() const { return accumulated_time_Heff + accumulated_time_Energy; };
    std::pair<double, double> elapsed_time_Heff_and_Energy() const {
        return {accumulated_time_Heff, accumulated_time_Energy};
    }

    /// For wrapping only: pointer to h()
    virtual long int h_ptr(const State& state) const { return (long int)(new af::array(H_eff(state)))->get(); }

  protected:
    ///< Calculating the micromagnetic energy from the h field
    virtual double E(const State& state, const af::array& h) const = 0;

  private:
    double E(const State& state) const { return Energy(state, H_eff(state)); };

    virtual af::array h(const State& state) const = 0;
    mutable double accumulated_time_Heff{0.};
    mutable double accumulated_time_Energy{0.};
};

constexpr bool timing_is_on{true};

inline af::array FieldTerm::H_eff(const State& state) const {
    if (timing_is_on) {
        af::timer timer = af::timer::start();
        const auto result = h(state);
        accumulated_time_Heff += timer.stop();
        return result;
    } else {
        return h(state);
    }
}

inline double FieldTerm::Energy(const State& state, const af::array& h) const {
    if (timing_is_on) {
        af::timer timer = af::timer::start();
        const auto result = E(state, h);
        accumulated_time_Energy += timer.stop();
        return result;
    } else {
        return E(state, h);
    }
}

inline double FieldTerm::Energy(const State& state) const {
    if (timing_is_on) {
        af::timer timer = af::timer::start();
        const auto result = E(state, H_eff(state));
        accumulated_time_Energy += timer.stop();
        return result;
    } else {
        return E(state, H_eff(state));
    }
}

// Aliases used to initialize objects wich inherit from this class
using uptr_FieldTerm = std::unique_ptr<FieldTerm>;
using vec_uptr_FieldTerm = std::vector<uptr_FieldTerm>;

namespace fieldterm {

// template <typename T> std::ostream& print_elapsed_time(std::ostream& os, const T& fieldterms) {
template <typename T> void print_elapsed_time(const T& fieldterms, std::ostream& os = std::cout) {
    for (const auto& elem : fieldterms) {
        const auto i = &elem - &fieldterms[0];
        const auto [H_eff, Energy] = elem->elapsed_time_Heff_and_Energy();
        os << "[" << i << "] elapsed time [s]: H_eff: " << H_eff << ", Energy: " << Energy << std::endl;
    }
}

/// Calculate effective field by accumulating all h(state) terms in container.
/// Expects non-empty container with at least one element
template <typename T> af::array accumulate_heff(const T& fieldterms, const State& state) {
    return std::accumulate(std::begin(fieldterms) + 1, std::end(fieldterms), fieldterms[0]->H_eff(state),
                           [&state](const auto& sum, const auto& elem) { return sum + elem->H_eff(state); });
}

/// Accumulates Energy(state) of all fieldterms
template <typename T> double accumulate_E(const T& fieldterms, const State& state) {
    return std::accumulate(std::begin(fieldterms), std::end(fieldterms), double{},
                           [&state](const auto& sum, const auto& elem) { return sum + elem->Energy(state); });
}

template <typename T, typename U> std::pair<T, U> operator+(const std::pair<T, U>& l, const std::pair<T, U>& r) {
    return {l.first + r.first, l.second + r.second};
}

/// Calculates accumulated effective field via h(state) and energy E via E(state,h).
/// Expects non-empy container.
/// makes use of E(state, heff) overload, avoiding second heff calculation when compared to E(state)
template <typename T> std::pair<af::array, double> accumulate_heff_and_E(const T& fieldterms, const State& state) {
    return std::accumulate(std::begin(fieldterms) + 1, std::end(fieldterms), fieldterms[0]->h_and_E(state),
                           [&state](const auto& sum, const auto& elem) { return sum + elem->h_and_E(state); });
}

// create a unique_ptr (e.g. from ctor or copy-ctor)
template <typename T, class... Args> auto fieldterm_ptr(Args&&... args) {
    return std::unique_ptr<FieldTerm>(std::make_unique<T>(std::forward<Args>(args)...));
}

// create unique_ptr to a copy
template <typename T> std::unique_ptr<FieldTerm> cp_to_uptr(const T& t) { return std::make_unique<T>(t); }

// moves to a unique_ptr, do not use element afterwards
// requires std::move() at callsite
template <typename T> std::unique_ptr<FieldTerm> mv_to_uptr(T&& t) { return std::make_unique<T>(t); }

/// returns std::vector<std::unique_ptr<FieldTerm>> from args
/// args called with std::move() are moved, copied elsewise
template <typename... Args> auto to_vec(Args... args) {
    std::vector<std::unique_ptr<FieldTerm>> v;
    (v.push_back(std::unique_ptr<FieldTerm>(std::make_unique<Args>(args))), ...);
    return v;
}

// Always moves arguments, except when arg is const itself
// Uses mv-ctor on args or copy-ctor as fallback (when arg is const)
template <typename T> std::unique_ptr<FieldTerm> to_uptr(T t) { return std::make_unique<T>(t); }
template <typename... Args> auto mv_to_vec(Args&&... args) {
    std::vector<std::unique_ptr<FieldTerm>> v;
    (v.push_back(to_uptr(std::move(args))), ...);
    return v;
}

} // namespace fieldterm
} // namespace magnumafcpp
