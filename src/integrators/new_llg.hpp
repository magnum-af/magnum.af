#ifndef NEW_LLG_H
#define NEW_LLG_H
#include "arrayfire.h"
#include "../state.hpp"
#include "../func.hpp"
#include "../llg_terms/LLGTerm.hpp"
#include "adaptive_runge_kutta.hpp"
#include <memory>//shared_ptr


class NewLlg {
    public:
        NewLlg(std::string scheme = "RKF45");
        void step(State&);
        double E(const State&);
        const bool dissipation_term_only{false};
        std::vector<std::shared_ptr<LLGTerm> > Fieldterms;
    private:
        af::array fdmdt(const State& state);
        typedef array (NewLlg::*fptr)(const State& state);
	fptr GetFPointer(void)
        {
            return &NewLlg::fdmdt;
        }
        af::array fheff(const State& state);
        AdaptiveRungeKutta integrator;
         
};

#endif
