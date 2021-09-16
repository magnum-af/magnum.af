#pragma once
#include "adaptive_runge_kutta.hpp"
#include "arrayfire.h"
#include "field_terms/field_term.hpp"
#include "state.hpp"

namespace magnumaf {
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
    LLGIntegrator(double alpha, vec_uptr_FieldTerm llgterms, const std::string& scheme = "RKF45",
                  Controller controller = Controller(), bool dissipation_term_only = false);
    LLGIntegrator(double alpha, std::initializer_list<movable_il<uptr_FieldTerm>> llgterms,
                  const std::string& scheme = "RKF45", Controller controller = Controller(),
                  bool dissipation_term_only = false);
    double alpha{0}; //!< Unitless damping constant in the
                     //!< Landau-Lifshitz-Gilbert equation

    vec_uptr_FieldTerm llgterms;

    double E(const State&) const;

    void relax(State& state, double precision = 1e-10, unsigned iloop = 100, unsigned iwritecout = 1000,
               bool verbose = true);

    /// Integrate and produce dense output.
    /// @param [in, out] state State which is to be integrated
    /// @param [in] time_in_s Time to integrate in [s]
    /// @param [in] write_every_dt_in_s Timestap for which t, mx, my, mz is printed
    /// @param [out] os outputstream
    /// @param [in] verbose Verbose switch
    void integrate_dense(State& state, double time_in_s, double write_every_dt_in_s, std::ostream& os = std::cout,
                         bool verbose = true);
    // wrapping using filename to omit stream:
    void integrate_dense(State& state, double time_in_s, double write_every_dt_in_s, const std::string& filename,
                         bool verbose = true, bool append = false);
    double get_time_heff() const { return time_heff; }
    long int h_addr(const State& state) const;

  private:
    af::array fheff(const State& state) const;
    const bool dissipation_term_only;
    virtual af::array f(const State& state) const override;
    mutable double time_heff{0};
};

} // namespace magnumaf
