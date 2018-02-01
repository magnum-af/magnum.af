#ifndef STOCHASTIC_INTEGRATOR_H
#define STOCHASTIC_INTEGRATOR_H
#include "arrayfire.h"
#include <memory>
#include <chrono>
#include "LLGTerm.hpp"
#include "state.hpp"
#include "func.hpp"
using namespace af;

class Stochastic_Integrator {
    public:
        Stochastic_Integrator (State, std::vector<std::shared_ptr<LLGTerm> >, const double, std::string );
        std::vector<std::shared_ptr<LLGTerm> > Fieldterms;
        Param param;
        Mesh mesh;
        void step(State& state, const double dt);
        double cpu_time();

        //Getter functions
        unsigned long int get_time() const { return time;};
        unsigned long int get_calls() const { return calls ;};
        unsigned long int get_fdmdt_calls() const { return fdmdt_calls ;};
        unsigned long int get_stochfdmdt_calls() const { return stochfdmdt_calls ;};

    protected:
        template <class T> T Heun(const T& m, const double dt);
        template <class T> T SemiImplicitHeun(const T& m, const double dt);
        template <class T> T detRK4(const T& m, const double dt);

        unsigned long int fdmdt_calls{0};
        unsigned long int calls{0};
        unsigned long int stochfdmdt_calls{0};
        double time{0.};
    private:
        virtual array detfdmdt(const array& m)=0;
        virtual array stochfdmdt(const array& m, const array& h_th)=0;

        array m_prev;
        array h_th_prev;
        af::timer timer_stoch;
       
        int mode; //Integration mode

        af::randomEngine rand_engine;
};

#endif
