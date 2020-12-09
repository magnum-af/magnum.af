#pragma once
#include "adaptive_runge_kutta.hpp"
#include "arrayfire.h"
#include "field_terms/field_term.hpp"
#include "state.hpp"

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
    LLGIntegrator(double alpha, vec_uptr_FieldTerm llgterms, std::string scheme = "RKF45",
                  Controller controller = Controller(), bool dissipation_term_only = false);
    LLGIntegrator(double alpha, std::initializer_list<movable_il<uptr_FieldTerm>> llgterms,
                  std::string scheme = "RKF45", Controller controller = Controller(),
                  bool dissipation_term_only = false);
    double alpha{0}; //!< Unitless damping constant in the
                     //!< Landau-Lifshitz-Gilbert equation

    vec_uptr_FieldTerm llgterms;

    double E(const State&) const;

    void relax(State& state, double precision = 1e-10, unsigned iloop = 100, unsigned iwritecout = 1000,
               bool verbose = true);

    double get_time_heff() const { return time_heff; }
    long int h_addr(const State& state) const;

    af::array fheff(const State& state) const; // TODO move to priv again

  private:
    const bool dissipation_term_only;
    virtual af::array f(const State& state) const override;
    mutable double time_heff{0};
};

} // namespace magnumafcpp
