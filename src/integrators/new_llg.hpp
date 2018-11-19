#ifndef NEW_LLG_H
#define NEW_LLG_H
#include "arrayfire.h"
#include "../state.hpp"
#include "../func.hpp"
#include "../llg_terms/LLGTerm.hpp"
#include "adaptive_runge_kutta.hpp"
#include <memory>//shared_ptr

class NewLlg : public AdaptiveRungeKutta{
    public:
        NewLlg(std::string scheme = "RKF45", Controller controller = Controller(), bool dissipation_term_only = false);
        NewLlg(LlgTerms llgterms, std::string scheme = "RKF45", Controller controller = Controller(), bool dissipation_term_only = false);
        double E(const State&);
        LlgTerms llgterms;
        const bool dissipation_term_only;
        double get_time_heff(){return time_heff;}
        void relax(State& state, double precision = 1e-10, const int iloop = 100, const int iwritecout = 1000);
        long int get_fheff_addr(const State& state);
        af::array fheff_tmp;//TODO tempfix for wrapping, elaborate other solution
    private:
        af::array f(const State& state);
        af::array fheff(const State& state);
        double time_heff{0};
         
};

#endif
