#include "arrayfire.h"
#include "magnum_af.hpp"
#include "integrators/integrator.hpp"
#include "integrators/adaptive_runge_kutta.hpp"

using namespace magnumaf;


using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;

void calcm(State state, std::ostream& myfile);

af::array givem (const float t, const af::array& m){
    return t*sqrt(m);
}
int main(int argc, char** argv)
{
    cout.precision(12);
    AdaptiveRungeKutta rkf(&givem, "BS45", Controller(1e-15, 1e4, 1e-10, 1e-10));
    float t=0;
    array m = constant(1.0, 1, f32);
    for (int i=0; i<1e4; i++){
         rkf.step(m, t);
         float mcalc=afvalue(m);
         float manal=1./16. * pow(pow(t, 2)+4, 2);
         std::cout << "t= "<< t << ", m= " << mcalc << ", analy= " << manal << ", rel_err= " << (mcalc-manal)/(mcalc+manal) <<  std::endl;
    }
}

void calcm(State state, std::ostream& myfile){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m, 0)<< "\t" <<meani(state.m, 1)<< "\t" <<meani(state.m, 2)<< "\t" << std::endl;
}
