#pragma once
#include "field_terms/field_term.hpp"
#include <fstream>
#include <utility>

namespace magnumaf {

/// LBFGS Minimizer, adapted from bvec implementation, courtesy of Thomas Schrefl.
class LBFGS_Minimizer {
  public:
    LBFGS_Minimizer(double tolerance = 1e-6, size_t maxIter = 230, int verbose = 4);
    LBFGS_Minimizer(vec_uptr_FieldTerm llgterms, double tolerance_ = 1e-6, size_t maxIter_ = 230, int verbose = 4);

    double Minimize(State&) const;

    vec_uptr_FieldTerm fieldterms_{}; // default init, as not constructed in init list

    mutable std::ofstream of_convergence_; // stream to write additional convergence data to

  private:
    const double tolerance_; ///< Error tolerance with default 1e-6
    const size_t maxIter_;   ///< Maximum number of iterations
    const int verbose_;      ///< Setting output options, valid values are 0, 1, 2, 3, 4
};

} // namespace magnumaf
