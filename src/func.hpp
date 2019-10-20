#pragma once
#include <iostream>
#include <chrono> //for af::randomEngine
#include "arrayfire.h"

namespace magnumaf{


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

af::array cross4(const af::array& a, const af::array& b);
af::array dotproduct(const af::array& a, const af::array& b);
af::array renormalize(const af::array& a);
af::array renormalize_handle_zero_values(const af::array& a);
af::array vecnorm(const af::array& a);
float afvalue(const af::array& a); //give value of a 1, 1, 1, 1 af af::array
unsigned int afvalue_u32(const af::array& a); // Returns value an af::array of type u32 == 6 and size [1, 1, 1, 1]
float full_inner_product(const af::array& a, const af::array& b);
float maxnorm(const af::array& a);
float minval(const af::array& a);
float meani(const af::array& a, const int i);
float FrobeniusNorm(const af::array& a);
//TODO void calcm(State state, LLG Llg, std::ostream& myfile);
float euclnorm(const af::array& a);
//TODO auto rk4(af::array f(float, af::array));
float mean_abs_diff(const af::array& a, const af::array& b); //!< absolute difference
float mean_rel_diff(const af::array& a, const af::array& b); //!< relative difference
float max_abs_diff(const af::array& a, const af::array& b); //!< absolute difference
float max_rel_diff(const af::array& a, const af::array& b); //!< relative difference
bool abs_diff_lt_precision(af::array first, af::array second, float precision = 4e-8, bool verbose = true); //!< Absolute difference less than precision
bool rel_diff_lt_precision(af::array first, af::array second, float precision = 2e-3, bool verbose = true); //!< Relative difference less than precision
float abs_diff_upperbound(const af::array& a, const af::array& b, bool verbose = true, float start_precision = 1e0, float factor1 = 0.1, float factor2 = 0.9);//!< Calculates upper bound with abs_diff_lt_precision
float rel_diff_upperbound(const af::array& a, const af::array& b, bool verbose = true, float start_precision = 1e0, float factor1 = 0.1, float factor2 = 0.9);//!< Calculates upper bound with rel_diff_lt_precision

// utility functions
namespace util{
    std::pair<int, int> k2ij(const int k, const int n);//!< returns indices i, j from a serialized triangular matrix index k
    int ij2k(const int i, const int j, const int n);//!< returns the serialized triangular matrix index k from regular indices i, j
    unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int ni);
    unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int ni, const unsigned int nj);
    unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l, const unsigned int ni, const unsigned int nj, const unsigned int nk);

    af::randomEngine rand_engine_current_time();
}
}// namespace magnumaf
