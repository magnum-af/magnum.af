#ifndef NEW_LLG_H
#define NEW_LLG_H
#include "arrayfire.h"
#include "../state.hpp"
#include "../func.hpp"
#include "../llg_terms/LLGTerm.hpp"
#include "adaptive_runge_kutta.hpp"
#include <memory>//shared_ptr

typedef std::shared_ptr<LLGTerm> LlgTerm; 
typedef std::vector<LlgTerm> LlgTerms; 

class NewLlg : public AdaptiveRungeKutta{
    public:
        NewLlg(std::string scheme = "RKF45", Controller controller = Controller(), bool dissipation_term_only = false);
        double E(const State&);
        const bool dissipation_term_only;
        LlgTerms llgterms;
        double get_time_heff(){return time_heff;}
    private:
        af::array f(const State& state);
        af::array fheff(const State& state);
        double time_heff{0};
         
};

#endif
