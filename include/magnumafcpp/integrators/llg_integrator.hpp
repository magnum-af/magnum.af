#pragma once
#include "field_terms/LLGTerm.hpp"
#include "adaptive_runge_kutta.hpp"
#include "arrayfire.h"

namespace magnumafcpp {
template <class T> struct movable_il {
    mutable T t;
    operator T() const&& { return std::move(t); }
    movable_il(T&& in) : t(std::move(in)) {}
};

///
/// The LLGIntegrator class performs time integration of the
/// Landau–Lifshitz–Gilbert (LLG) equation on the magnetization passed by an
/// State object. The LLG equation reads \f[
///     \frac{\partial \boldsymbol{m}}{\partial t} =
///             - \gamma \boldsymbol{m} \times \boldsymbol{H}_{\text{eff}}
///             - \alpha \gamma \boldsymbol{m} \times (\boldsymbol{m} \times
///             \boldsymbol{H}_{\text{eff}})
/// \f]
///

class LLGIntegrator : public AdaptiveRungeKutta {
  public:
    LLGIntegrator(double alpha, std::string scheme = "RKF45", Controller controller = Controller(),
                  bool dissipation_term_only = false);
    LLGIntegrator(double alpha, LlgTerms llgterms, std::string scheme = "RKF45", Controller controller = Controller(),
                  bool dissipation_term_only = false);
    LLGIntegrator(double alpha, std::initializer_list<movable_il<LlgTerm>> llgterms, std::string scheme = "RKF45",
                  Controller controller = Controller(), bool dissipation_term_only = false);
    double alpha{0}; //!< Unitless damping constant in the
                     //!< Landau-Lifshitz-Gilbert equation
    LlgTerms llgterms;
    const bool dissipation_term_only;
    double E(const State&);

    double get_time_heff() { return time_heff; }
    void relax(State& state, double precision = 1e-10, const unsigned iloop = 100, const unsigned iwritecout = 1000,
               const bool verbose = true);
    long int h_addr(const State& state);

  private:
    af::array f(const State& state);
    af::array fheff(const State& state);
    double time_heff{0};
};

} // namespace magnumafcpp
