#ifndef NEW_LLG_H
#define NEW_LLG_H
#include "arrayfire.h"
#include "../state.hpp"
#include "../llg_terms/LLGTerm.hpp"


class NewLlg {
   public:
        void step(State&);
        double E(const State&);
        const dissipation_term_only{false};
   private:
         
};

#endif
