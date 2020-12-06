#pragma once
#include "arrayfire.h"
#include "field_terms/field_term.hpp"
#include "state.hpp"
#include <memory>

namespace magnumafcpp {

class Stochastic_Integrator {
  public:
    Stochastic_Integrator(double alpha, double T, double dt, State state,
                          std::vector<std::unique_ptr<FieldTerm>> fieldterms_in, std::string smode);
    virtual ~Stochastic_Integrator() = default;
    double alpha;
    double T;  //<! Temparature in [K]
    double dt; //<! Timestep in [s]
    std::vector<std::unique_ptr<FieldTerm>> fieldterms;
    void step(State& state);
    double cpu_time();

    // Getter functions
    unsigned long int get_calls() const { return calls; };
    unsigned long int get_fdmdt_calls() const { return fdmdt_calls; };
    unsigned long int get_stochfdmdt_calls() const { return stochfdmdt_calls; };
    double get_time() const { return timer; };

  protected:
    af::array Heun(const State&);
    af::array SemiImplicitHeun(const State&);
    af::array detRK4(const State&);

    unsigned long int calls{0};
    unsigned long int fdmdt_calls{0};
    unsigned long int stochfdmdt_calls{0};
    double timer{0.};

  private:
    virtual af::array detfdmdt(const State&) = 0;
    virtual af::array stochfdmdt(const State&, const af::array& h_th) = 0;

    af::array m_prev;
    af::array h_th_prev{};

    const int mode; // Integration mode

    af::randomEngine rand_engine;
};

} // namespace magnumafcpp
