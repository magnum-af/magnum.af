#pragma once
#include "../state.hpp"
#include "../constants.hpp"
#include "../func.hpp"
#include "arrayfire.h"
#include <memory>

namespace magnumafcpp{

// Abstract basis class for all terms in the LLG equation.
class LLGTerm{
  public:
    virtual af::array h (const State& state) =0;
    /// Calculating the micromagnetic energy \f$E\f$.
    /// This is a prototype for all llgterms with are linear in m and must be overwritten in e.g. zeeman where factor 1/2 becomes 1.
    virtual double E (const State& state){
        if( state.Ms_field.isempty() ){
            return -constants::mu0/2. * state.Ms * afvalue(af::sum(af::sum(af::sum(af::sum( h(state) * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
        }
        else{
            return -constants::mu0/2. * afvalue(af::sum(af::sum(af::sum(af::sum(state.Ms_field * h(state) * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
        }
    }
    ///< Calculating the micromagnetic energy for a already calculated h field (to save computational cost)
    virtual double E (const State& state, const af::array& h){
        if( state.Ms_field.isempty() ){
            return -constants::mu0/2. * state.Ms * afvalue(af::sum(af::sum(af::sum(af::sum( h * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
        }
        else{
            return -constants::mu0/2. * afvalue(af::sum(af::sum(af::sum(af::sum(state.Ms_field * h * state.m, 0), 1), 2), 3)) * state.mesh.dx * state.mesh.dy * state.mesh.dz;
        }
    }
    virtual double get_cpu_time()=0;

    /// For wrapping only: pointer to h()
    virtual long int h_ptr(const State& state){
        return (long int) (new af::array(h(state)))->get();
    }
    virtual ~LLGTerm(){};
};


// Aliases used to initialize objects wich inherit from this class
using LlgTerm = std::shared_ptr<LLGTerm>;
using LlgTerms = std::vector<LlgTerm>;

}// namespace magnumafcpp
