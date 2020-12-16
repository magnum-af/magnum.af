#pragma once
#include "arrayfire.h"
#include <array>
#include <chrono> //for af::randomEngine
#include <iostream>

namespace magnumafcpp {

class WrappedArray {
  public:
    WrappedArray(af::array array);
    WrappedArray(long int array_ptr);
    ~WrappedArray(){};
    af::array array{};
    void set_array(long int array_ptr);
    long int get_array_addr();

  private:
};

///
/// Function calculating the spacial mean (i.e. mean along first three
/// dimensions) in a specified region only of a given vectorfield with size [nx,
/// ny, nz, 3]. The region is specified by an array of size [nx, ny, nz, 1],
/// only cells with non-zero values are considered for the mean. expects sizes
/// to be:
///     vectorfield [nx, ny, nz, 3]
///     region      [nx, ny, nz, 1]
///
std::array<double, 3> spacial_mean_in_region(const af::array& vectorfield, const af::array& region);
std::array<double, 3> spacial_mean_in_region(long int vectorfield,
                                             long int region); // Wrapping only

af::array cross4(const af::array& a, const af::array& b);
af::array cross4shift(const af::array& a, const af::array& b);
af::array dotproduct(const af::array& a, const af::array& b);
af::array normalize(const af::array& a);
af::array normalize_handle_zero_vectors(const af::array& a);
void normalize_inplace(af::array& a);
af::array vecnorm(const af::array& a);
double afvalue(const af::array& a);           // give value of a 1, 1, 1, 1 af af::array
unsigned int afvalue_u32(const af::array& a); // Returns value an af::array of type u32 == 6
                                              // and size [1, 1, 1, 1]
double full_inner_product(const af::array& a, const af::array& b);
double min_4d(const af::array& a);
double meani(const af::array& a, const int i);
double FrobeniusNorm(const af::array& a);
// TODO void calcm(State state, LLG llg, std::ostream& myfile);
double euclnorm(const af::array& a);
// TODO auto rk4(af::array f(double, af::array));
double mean_abs_diff(const af::array& a,
                     const af::array& b); //!< absolute difference
double mean_rel_diff(const af::array& a,
                     const af::array& b); //!< relative difference
double max_abs_diff(const af::array& a,
                    const af::array& b); //!< absolute difference
double max_rel_diff(const af::array& a,
                    const af::array& b); //!< relative difference
bool abs_diff_lt_precision(af::array first, af::array second, double precision = 4e-8,
                           bool verbose = true); //!< Absolute difference less than precision
bool rel_diff_lt_precision(af::array first, af::array second, double precision = 2e-3,
                           bool verbose = true); //!< Relative difference less than precision
double abs_diff_upperbound(const af::array& a, const af::array& b, bool verbose = true, double start_precision = 1e0,
                           double factor1 = 0.1,
                           double factor2 = 0.9); //!< Calculates upper bound with abs_diff_lt_precision
double rel_diff_upperbound(const af::array& a, const af::array& b, bool verbose = true, double start_precision = 1e0,
                           double factor1 = 0.1,
                           double factor2 = 0.9); //!< Calculates upper bound with rel_diff_lt_precision

// utility functions
namespace util {
std::pair<int, int> k2ij(const int k,
                         const int n); //!< returns indices i, j from a
                                       //!< serialized triangular matrix index k
int ij2k(const int i, const int j,
         const int n); //!< returns the serialized triangular matrix index k
                       //!< from regular indices i, j
unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int ni);
unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int ni,
                    const unsigned int nj);

/// Get the flat index from af::array with 4d indices, i.e. af::array(i, j, k, l), where af::dims=[ni, nj, nk, :].
unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l,
                    const unsigned int ni, const unsigned int nj, const unsigned int nk);

af::randomEngine rand_engine_current_time();
} // namespace util
} // namespace magnumafcpp
