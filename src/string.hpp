#pragma once
#include "state.hpp"
#include "integrators/new_llg.hpp"
#include <vector>
#include <algorithm>

namespace magnumaf{


class String {
  public:
    String(float alpha, State state, std::vector<State> inputimages, int n_interp, float dt, LlgTerms Fieldterms_in);
    float run(const std::string filepath, const float string_abort_rel_diff = 1e-12, const float string_abort_abs_diff = 1e-27, const int string_steps = 10000);

    State state;// TODO remove this state instace
    //Mesh mesh;
    LLGIntegrator Llg;//(state_relax, atol, rtol, hmax, hmin);
    const int n_interp;
    const float dt;
    float time{0};
    std::vector<float> x;//Current x values
    std::vector<float> x_interp;//x-values where to interpolate at (will be regular grid)
    std::vector<float> E;//Energy values
    std::vector<State> images;//Current images
    void calc_E();
    void calc_x();
    void calc_x(std::vector<State>);
    void lin_interpolate();
    void integrate();//Integrate all images for dt
    void step();
    void vec_renormalize();
    void write_vti(std::string file);
  private:
};

}// namespace magnumaf
