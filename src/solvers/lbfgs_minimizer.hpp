#pragma once
#include "../llg_terms/LLGTerm.hpp"
#include <fstream>

namespace magnumaf{


//For second Method, use interface class:
//https://stackoverflow.com/questions/40624175/c-how-to-implement-a-switch-between-class-members
//
//typedef float Dtype;// TODO: replace ?
  //NOTE: Dtype ak in Schrefl::linesearch is moved into cvsrch

class LBFGS_Minimizer {
    public:
        LBFGS_Minimizer(float tolerance_ = 1e-6, size_t maxIter_ = 230, int verbose = 4);
        LBFGS_Minimizer(LlgTerms llgterms, float tolerance_ = 1e-6, size_t maxIter_ = 230, int verbose = 4);

        float Minimize(State&);

        LlgTerms llgterms_;

        float GetTimeCalcHeff() const { return time_calc_heff_;}; ///< Accumulated time for calculation of Heff.
        std::ofstream of_convergence;
    private:
        af::array Gradient(const State&);///< Calculate gradient as energy-dissipation term of llg
        af::array Heff(const State& m);///< Effective Field
        float E(const State&); ///< Calculate Energy
        float EnergyAndGradient(const State& state, af::array& gradient);
        float time_calc_heff_{0};///< Timer measuring calls to effective field _h
        const float tolerance_;///< Error tolerance with default 1e-6
        const size_t maxIter_;///< Maximum number of iterations
        const int verbose;///< Setting output options, valid values are 0, 1, 2, 3, 4
        float mxmxhMax(const State& state);///< TODO investigate definition, init value etc
        float linesearch(State& state, float &fval, const af::array &x_old, af::array &g, const af::array &searchDir, const float tolf);
        //TODO//TODEL//int cvsrch(const State& state, const af::array &wa, af::array &x, float &f, af::array &g, const af::array &s, float tolf);
        int cvsrch(State& state, const af::array &wa, float &f, af::array &g, float &stp, const af::array &s, const float tolf);
        int cstep(float& stx, float& fx, float& dx, float& sty, float& fy, float& dy, float& stp, float& fp, float& dp, bool& brackt, float& stpmin, float& stpmax, int& info);
};

}// namespace magnumaf
