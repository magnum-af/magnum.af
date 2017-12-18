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
    template <class T>  T test(const T& m, const double dt);
    Stochastic_LLG (State, std::vector<std::shared_ptr<LLGTerm> >);
    std::vector<std::shared_ptr<LLGTerm> > Fieldterms;

    //array rk4(const array& m, const double dt);
    template <class T>  T rk4(const T& m, const double dt);
    array SemiHeun(const array& m, const double dt);
    template <class T> T SemiImplicitHeun(const T& m, const double dt);
    template <class T> T StochSemiImplicitHeun(const T& m, const double dt);
    //void SemiImplicitHeun(array& m, const double dt);
    //array SemiImplicitHeun(const array& m, const double dt);
    void step(State& state, const double dt);
    unsigned long int get_calls() const { return calls ;};
    unsigned long int get_fdmdt_calls() const { return fdmdt_calls ;};
    unsigned long int get_stochfdmdt_calls() const { return stochfdmdt_calls ;};
    double cpu_time();

    Param param;
    double     time{0.};
    private:
    unsigned long int calls{0};
    unsigned long int fdmdt_calls{0};
    unsigned long int stochfdmdt_calls{0};
    Mesh mesh;
    array fheff(const array& m);
    array fdmdt(const array& m);
    array stochfdmdt(const array& m, const array& h_th);

    array m_prev;
    af::timer timer_stoch;
};

#endif
