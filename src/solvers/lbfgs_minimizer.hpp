#ifndef LBFGS_MINIMIZER_H
#define LBFGS_MINIMIZER_H
#include <memory>
#include <list>
#include <algorithm>
#include "arrayfire.h"
#include "../state.hpp"
#include "../misc.hpp"
#include "../func.hpp"
#include "../llg_terms/LLGTerm.hpp"

//For second Method, use interface class:
//https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members
//
//typedef double Dtype;// TODO: replace ?
  //NOTE: Dtype ak in Schrefl::linesearch is moved into cvsrch

class LBFGS_Minimizer {
    public:
        LBFGS_Minimizer();

        double Minimize(State&);

        LlgTerms llgterms_;

        double GetTimeCalcHeff() const { return time_calc_heff_;}; ///< Accumulated time for calculation of Heff.
    private:
        af::array Gradient(const State&);///< Calculate gradient as energy-dissipation term of llg
        af::array Heff(const State& m);///< Effective Field 
        double E(const State&); ///< Calculate Energy
        double time_calc_heff_{0};///< Timer measuring calls to effective field _h
        int verbose_{5};///< Setting output options
        size_t maxIter_{200};///< Maximum number of iterations
        double mxmxhMax(const State& state);///< TODO investigate definition, init value etc
        double linesearch(const State& state, double &fval, const af::array &x_old, af::array &x, af::array &g, const af::array &searchDir, double tolf);
        //TODO//TODEL//int cvsrch(const State& state, const af::array &wa, af::array &x, double &f, af::array &g, const af::array &s, double tolf);
        int cvsrch(const State& state, const af::array &wa, af::array &x, double &f, af::array &g, double &stp, const af::array &s, double tolf);
        int cstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, double& stp, double& fp, double& dp, bool& brackt, double& stpmin, double& stpmax, int& info);
};

#endif
