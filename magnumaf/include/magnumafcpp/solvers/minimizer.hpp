#pragma once
#include "arrayfire.h"
#include "field_terms/field_term.hpp"
#include "state.hpp"

namespace magnumaf {

// For second Method, use interface class:
// https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members

// Energy minimizer using a semi-implicit update scheme applying the
// Barzilian-Borwein (BB) rule for stepsize calculation.

class Minimizer {
  public:
    Minimizer(std::string scheme = "BB", double tau_min = 1e-10, double tau_max = 1e-5, double dm_max = 1e4,
              int samples = 10, bool info = false);

    void minimize(State&) const; // Minimization routine

  private:
    af::array h(const State& m) const; // Effective Field
    af::array dm(const State& state) const;
    af::array m_next(const State& state, const double tau) const;
    double E(const State& state) const; // only for testing, remove?

    vec_uptr_FieldTerm fieldterms;

    std::string scheme;
    double tau_min;
    double tau_max;
    double dm_max;
    unsigned int samples;
    bool info; // if ture, prints step info to cout
};

} // namespace magnumaf
