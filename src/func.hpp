#ifndef FUNC_H
#define FUNC_H
#include <iostream>
#include "arrayfire.h"

class WrappedArray {
    public:
        WrappedArray(af::array array);
        WrappedArray(long int array_ptr);
        ~WrappedArray(){};
        af::array array;
        void set_array(long int array_ptr);
        long int get_array_addr();
    private:
};

af::array cross4(const af::array& a,const af::array& b);
af::array dotproduct(const af::array& a,const af::array& b);
af::array renormalize(const af::array& a);
af::array renormalize_handle_zero_values(const af::array& a);
af::array vecnorm(const af::array& a);
double afvalue(const af::array& a); //give value of a 1,1,1,1 af af::array
unsigned int afvalue_u32(const af::array& a); // Returns value an af::array of type u32 == 6 and size [1,1,1,1]
double full_inner_product(const af::array& a, const af::array& b);
double maxnorm(const af::array& a);
double minval(const af::array& a);
double meani(const af::array& a, const int i);
double FrobeniusNorm(const af::array& a);
//TODO void calcm(State state, LLG Llg, std::ostream& myfile);
double euclnorm(const af::array& a);
//TODO auto rk4(af::array f(double, af::array));
bool abs_diff_lt_precision(af::array first, af::array second, double precision = 4e-8, bool verbose = true); //!< Absolute difference less than precision
bool rel_diff_lt_precision(af::array first, af::array second, double precision = 2e-3, bool verbose = true); //!< Relative difference less than precision
double abs_diff_upperbound(const af::array& a, const af::array& b, bool verbose = true, double start_precision = 1e0, double factor1 = 0.1, double factor2 = 0.9);
double rel_diff_upperbound(const af::array& a, const af::array& b, bool verbose = true, double start_precision = 1e0, double factor1 = 0.1, double factor2 = 0.9);
#endif
