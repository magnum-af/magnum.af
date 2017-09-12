#ifndef STRING_H
#define STRING_H
#include <vector>
#include "arrayfire.h"
#include "llg.hpp"
class String {
  public:
    String(State state, std::vector<State> inputimages, int n_interp, double dt, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in);

    State state;
    //Mesh mesh;
    LLG Llg;//(state_relax,atol,rtol,hmax,hmin);
    int n_interp;
    double dt;
    double time{0};
    std::vector<double> x;//Current x values
    std::vector<double> x_interp;//x-values where to interpolate at (will be regular grid)
    std::vector<double> E;//Energy values
    std::vector<State> images;//Current images
    void calc_E();
    void calc_x();
    void calc_x(std::vector<State>);
    void lin_interpolate();
    void integrate();//Integrate all images for dt
    void step();
    void vec_renormalize();
  private:
};
#endif
