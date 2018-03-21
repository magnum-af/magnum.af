#ifndef FUNC_H
#define FUNC_H
#include "arrayfire.h"

////for calcm:
//#include "state.hpp"
//#include "llg.hpp"
//#include <iomanip>

using namespace af;
array cross4(const array& a,const array& b);
array dotproduct(const array& a,const array& b);
array renormalize(const array& a);
array renormalize_handle_zero_values(const array& a);
array vecnorm(const array& a);
double afvalue(const array& a); //give value of a 1,1,1,1 af array
double maxnorm(const array& a);
double meani(const array& a, const int i);
double FrobeniusNorm(const array& a);
//TODO void calcm(State state, LLG Llg, std::ostream& myfile);
double euclnorm(const array& a);
//TODO auto rk4(array f(double, array));
#endif
