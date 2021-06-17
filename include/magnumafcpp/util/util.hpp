#pragma once
#include "arrayfire.h"
#include <algorithm>
#include <numeric>

namespace magnumafcpp::util {

namespace pywrap {

/// Copies an af::array from python.
/// Make an af::array from a pointer obtained from python
inline af::array make_copy_form_py(long int py_ptr_to_afarray) {
    return *(new af::array(*((void**)py_ptr_to_afarray)));
}

/// Passes a copy of the passed array to python.
/// Allocates a copy of passed af::array, returns the af_array handle as a long int
// alternative technical name: get_ptr_to_copy
inline long int send_copy_to_py(const af::array& afarray) { return (long int)(new af::array(afarray))->get(); }

struct WrappedArray {
    explicit WrappedArray(af::array array) : array(std::move(array)) {}
    explicit WrappedArray(long int array_ptr) : array(make_copy_form_py(array_ptr)) {}
    af::array array{};
    void set_array(long int array_ptr) { this->array = make_copy_form_py(array_ptr); }
    long int get_array_copy_as_ptr() const { return send_copy_to_py(this->array); }
};

} // namespace pywrap

/// template function calculating the mean and the standard deviation of a given
/// template container data standard deviation with -1
// from
// https://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
template <class T> std::pair<double, double> mean_stdev_w_minus(T vec) {
    double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);
    double m = sum / vec.size();
    double accum = 0.0;
    std::for_each(std::begin(vec), std::end(vec), [&](const double d) { accum += (d - m) * (d - m); });
    double stdev = std::sqrt(accum / (vec.size() - 1));
    return {m, stdev};
}

/// standard deviation without -1
template <class T> std::pair<double, double> mean_stdev_no_minus(T vec) {
    double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);
    double m = sum / vec.size();
    double accum = 0.0;
    std::for_each(std::begin(vec), std::end(vec), [&](const double d) { accum += (d - m) * (d - m); });
    double stdev = std::sqrt(accum / (vec.size()));
    return {m, stdev};
}

inline std::array<double, 3> cross_product(std::array<double, 3> a, std::array<double, 3> b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

inline double dot_product(std::array<double, 3> a, std::array<double, 3> b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline double vector_norm(std::array<double, 3> vector) {
    return std::sqrt(std::pow(vector[0], 2) + std::pow(vector[1], 2) + std::pow(vector[2], 2));
}

inline std::array<double, 3> normalize_vector(std::array<double, 3> vector) {
    double norm = util::vector_norm(vector);
    return {vector[0] / norm, vector[1] / norm, vector[2] / norm};
}

template <typename T> std::array<T, 3> get_vec_scalars(const af::array& afvec) {
    const T xval = afvec(0, 0, 0, 0).scalar<T>();
    const T yval = afvec(0, 0, 0, 1).scalar<T>();
    const T zval = afvec(0, 0, 0, 2).scalar<T>();
    return {xval, yval, zval};
}

///
/// Function calculating the spacial mean (i.e. mean along first three
/// dimensions) in a specified region only of a given vectorfield with size [nx,
/// ny, nz, 3]. The region is specified by an array of size [nx, ny, nz, 1],
/// only cells with non-zero values are considered for the mean. expects sizes
/// to be:
///     vectorfield [nx, ny, nz, 3]
///     region      [nx, ny, nz, 1]
/// returns:
///     mean [1, 1, 1, 3]
///
af::array spacial_mean_in_region_afarray(const af::array& vectorfield, const af::array& region);

inline std::array<double, 3> spacial_mean_in_region(const af::array& vectorfield, const af::array& region) {
    return get_vec_scalars<double>(spacial_mean_in_region_afarray(vectorfield, region).as(f64));
}

// Wrapping only
inline std::array<double, 3> spacial_mean_in_region(long int vectorfield, long int region) {
    return spacial_mean_in_region(pywrap::make_copy_form_py(vectorfield), pywrap::make_copy_form_py(region));
}

af::array normalize(const af::array& a);
af::array normalize_handle_zero_vectors(const af::array& a);
void normalize_inplace(af::array& a);
af::array vecnorm(const af::array& a);
double afvalue_as_f64(const af::array& a); // give value of a 1, 1, 1, 1 af af::array
// Returns value an af::array of type u32 == 6
// and size [1, 1, 1, 1]
inline unsigned int afvalue_u32(const af::array& a) { return a.scalar<unsigned int>(); }
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

} // namespace magnumafcpp::util
