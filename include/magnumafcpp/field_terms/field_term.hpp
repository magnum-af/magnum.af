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
    af::array H_in_Apm(const State& state) const; ///< Effective field in Ampere per meter [A/m]
    af::array H_in_T(const State& state) const {
        return conversion::Apm_to_Tesla(H_in_Apm(state));
    }; ///< Effective field in Tesla [T]
    double Energy_in_J(const State& state, const af::array& h) const;
    double Energy_in_J(const State& state) const;
    double Energy_in_eV(const State& state, const af::array& h) const {
        return conversion::J_to_eV(Energy_in_J(state, h));
    };
    double Energy_in_eV(const State& state) const { return conversion::J_to_eV(Energy_in_J(state)); };

    /// Calculating the micromagnetic energy \f$E\f$.
    // virtual double impl_E_in_J(const State& state) const = 0;

    std::pair<af::array, double> H_in_Apm_Energy_in_J(const State& state) const {
        const auto htmp = H_in_Apm(state);
        return {htmp, Energy_in_J(state, htmp)};
    };

    double elapsed_eval_time() const { return accumulated_time_Heff + accumulated_time_Energy; };

    std::pair<double, double> elapsed_eval_time_Heff_and_Energy() const {
        return {accumulated_time_Heff, accumulated_time_Energy};
    }

    /// For wrapping only: raw pointer to copy of H_in_Apm(state)
    long int _pywrap_H_in_Apm(const State& state) const { return (long int)(new af::array(H_in_Apm(state)))->get(); }

  protected:
    ///< Calculating the micromagnetic energy from the h field
    virtual double impl_E_in_J(const State& state, const af::array& h) const = 0;

  private:
    double impl_E_in_J(const State& state) const { return Energy_in_J(state, H_in_Apm(state)); };

    virtual af::array impl_H_in_Apm(const State& state) const = 0;
    mutable double accumulated_time_Heff{0.};
    mutable double accumulated_time_Energy{0.};
};

constexpr bool timing_is_on{true};

inline af::array FieldTerm::H_in_Apm(const State& state) const {
    if (timing_is_on) {
        af::timer timer = af::timer::start();
        const auto result = impl_H_in_Apm(state);
        accumulated_time_Heff += timer.stop();
        return result;
    } else {
        return impl_H_in_Apm(state);
    }
}

inline double FieldTerm::Energy_in_J(const State& state, const af::array& h) const {
    if (timing_is_on) {
        af::timer timer = af::timer::start();
        const auto result = impl_E_in_J(state, h);
        accumulated_time_Energy += timer.stop();
        return result;
    } else {
        return impl_E_in_J(state, h);
    }
}

inline double FieldTerm::Energy_in_J(const State& state) const {
    if (timing_is_on) {
        af::timer timer = af::timer::start();
        const auto result = impl_E_in_J(state, H_in_Apm(state));
        accumulated_time_Energy += timer.stop();
        return result;
    } else {
        return impl_E_in_J(state, H_in_Apm(state));
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
        const auto [H_in_Apm, Energy] = elem->elapsed_time_Heff_and_Energy();
        os << "[" << i << "] elapsed time [s]: H_in_Apm: " << H_in_Apm << ", Energy: " << Energy << std::endl;
    }
}

/// Calculate effective field by accumulating all h(state) terms in container.
/// Expects non-empty container with at least one element
template <typename T> af::array Heff_in_Apm(const T& fieldterms, const State& state) {
    return std::accumulate(std::begin(fieldterms) + 1, std::end(fieldterms), fieldterms[0]->H_in_Apm(state),
                           [&state](const auto& sum, const auto& elem) { return sum + elem->H_in_Apm(state); });
}
template <typename T> af::array Heff_in_T(const T& fieldterms, const State& state) {
    return conversion::Apm_to_Tesla(Heff_in_Apm(fieldterms, state));
}

/// Accumulates Energy_in_J(state) of all fieldterms
template <typename T> double Eeff_in_J(const T& fieldterms, const State& state) {
    return std::accumulate(std::begin(fieldterms), std::end(fieldterms), double{},
                           [&state](const auto& sum, const auto& elem) { return sum + elem->Energy_in_J(state); });
}

template <typename T> double Eeff_in_eV(const T& fieldterms, const State& state) {
    return conversion::J_to_eV(Eeff_in_J(fieldterms, state));
}

template <typename T, typename U> std::pair<T, U> operator+(const std::pair<T, U>& l, const std::pair<T, U>& r) {
    return {l.first + r.first, l.second + r.second};
}

/// Calculates accumulated effective field via h(state) and energy E via E(state,h).
/// Expects non-empy container.
/// makes use of E(state, heff) overload, avoiding second heff calculation when compared to E(state)
template <typename T> std::pair<af::array, double> Heff_in_Apm_and_E(const T& fieldterms, const State& state) {
    return std::accumulate(
        std::begin(fieldterms) + 1, std::end(fieldterms), fieldterms[0]->H_in_Apm_Energy_in_J(state),
        [&state](const auto& sum, const auto& elem) { return sum + elem->H_in_Apm_Energy_in_J(state); });
}

// create a unique_ptr (e.g. from ctor or copy-ctor)
template <typename T, class... Args> auto to_uptr(Args&&... args) {
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
