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

    // void relax_dense(State& state, const vec_uptr_FieldTerm & fieldterms, size_t print_every_ns = 10){
    //}
    void run(double time_in_s, double print_every_s, State& state, const vec_uptr_FieldTerm& fieldterms,
             std::ostream& os = std::cout) {
        const double time_init = state.t;
        os << state << '\n';

        double time_last_write = state.t;
        double time_next_write = state.t;
        while (state.t < time_init + time_in_s) {
            step(state);
            if (state.t >= time_next_write) {
                // auto m_int = interpolate_m_at(time_next_write);
                // os << m_int << '\n';
                time_last_write = time_next_write;

                do {
                    time_next_write += print_every_s;
                } while (time_next_write < state.t); // in case state.t is already bigger, we skipp timestap
            }
        }
    }
    double get_time_heff() const { return time_heff; }
    long int h_addr(const State& state) const;

    af::array fheff(const State& state) const; // TODO move to priv again

  private:
    const bool dissipation_term_only;
    virtual af::array f(const State& state) const override;
    mutable double time_heff{0};
};

} // namespace magnumafcpp
