#pragma once
#include "../llg_terms/LLGTerm.hpp"
#include "adaptive_runge_kutta.hpp"
#include "arrayfire.h"

namespace magnumaf{

///
/// The LLGIntegrator class performs time integration of the Landau–Lifshitz–Gilbert (LLG) equation on the magnetization passed by an State object.
/// The LLG equation reads
/// \f[
///     \frac{\partial \boldsymbol{m}}{\partial t} =
///             - \gamma \boldsymbol{m} \times \boldsymbol{H}_{\text{eff}}
///             - \alpha \gamma \boldsymbol{m} \times (\boldsymbol{m} \times \boldsymbol{H}_{\text{eff}})
/// \f]
///

class LLGIntegrator : public AdaptiveRungeKutta{
    public:
        LLGIntegrator(float alpha, std::string scheme = "RKF45", Controller controller = Controller(), bool dissipation_term_only = false);
        LLGIntegrator(float alpha, LlgTerms llgterms, std::string scheme = "RKF45", Controller controller = Controller(), bool dissipation_term_only = false);
        float alpha{0};//!< Unitless damping constant in the Landau-Lifshitz-Gilbert equation
        LlgTerms llgterms;
        const bool dissipation_term_only;
        float E(const State&);

        float get_time_heff(){return time_heff;}
        void relax(State& state, float precision = 1e-10, const int iloop = 100, const int iwritecout = 1000);
        long int h_addr(const State& state);
    private:
        af::array f(const State& state);
        af::array fheff(const State& state);
        float time_heff{0};

};

}// namespace magnumaf
