#ifndef FUNC_H
#define FUNC_H
#include "arrayfire.h"

af::array cross4(const af::array& a,const af::array& b);
af::array dotproduct(const af::array& a,const af::array& b);
af::array renormalize(const af::array& a);
af::array renormalize_handle_zero_values(const af::array& a);
af::array vecnorm(const af::array& a);
double afvalue(const af::array& a); //give value of a 1,1,1,1 af af::array
double full_inner_product(const af::array& a, const af::array& b);
double maxnorm(const af::array& a);
double minval(const af::array& a);
double meani(const af::array& a, const int i);
double FrobeniusNorm(const af::array& a);
//TODO void calcm(State state, LLG Llg, std::ostream& myfile);
double euclnorm(const af::array& a);
//TODO auto rk4(af::array f(double, af::array));
#endif
