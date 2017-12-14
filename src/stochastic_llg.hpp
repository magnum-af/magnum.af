#ifndef STOCHASTIC_LLG_H
#define STOCHASTIC_LLG_H
#include "arrayfire.h"
#include <memory>
#include "LLGTerm.hpp"
#include "state.hpp"
#include "func.hpp"
using namespace af;

class Stochastic_LLG {
    public:
    Stochastic_LLG (State, std::vector<std::shared_ptr<LLGTerm> >);
    std::vector<std::shared_ptr<LLGTerm> > Fieldterms;

    array rk4(const array& m, const double dt);
    array StemiImplicitHeun(const array& m, const double dt);
    void step(State& state, const double dt);
    unsigned long int get_calls() const { return calls ;};
    unsigned long int get_fdmdt_calls() const { return fdmdt_calls ;};

    Param param;
    private:
    unsigned long int calls{0};
    unsigned long int fdmdt_calls{0};
    Mesh mesh;
    array fheff(const array& m);
    array fdmdt(const array& m);

    array m_prev;
};

#endif
