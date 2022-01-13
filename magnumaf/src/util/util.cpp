#include "util/util.hpp"
#include <chrono>
#include <iostream>

namespace magnumaf::util {

af::array spacial_mean_in_region_afarray(const af::array& vectorfield, const af::array& region) {
    const af::array zero_if_input_is_zero_else_one = !af::iszero(region);
    const auto number_of_nonzero_elements =
        af::sum(af::sum(af::sum(zero_if_input_is_zero_else_one, 0), 1), 2).scalar<unsigned>();
    const af::array considered_values = vectorfield * af::tile(zero_if_input_is_zero_else_one, 1, 1, 1, 3);
    const af::array sum_considered_values = af::sum(af::sum(af::sum(considered_values, 0), 1), 2);
    af::array mean_considered_values = sum_considered_values / number_of_nonzero_elements;
    return mean_considered_values;
}

/// Returns the value of array with only one element
double afvalue_as_f64(const af::array& a) {
    if (a.dims(0) != 1 || a.dims(1) != 1 || a.dims(2) != 1 || a.dims(3) != 1) {
        std::cout << "\33[1;31mWarning:\33[0m afvalue requested from array "
                     "with dim4 =["
                  << a.dims()
                  << "] != [1 1 1 1]. Only first entry will be returned. This "
                     "may lead to unexpected behaviour."
                  << std::endl;
    }
    return a.as(f64).scalar<double>();
}

double full_inner_product(const af::array& a, const af::array& b) {
    return afvalue_as_f64(sum(sum(sum(sum(a * b, 3), 2), 1), 0));
}

af::array normalize_even_zero_vectors(const af::array& a) { return a / tile(sqrt(sum(a * a, 3)), 1, 1, 1, 3); }

// Normalize all vectors to unit length, empty vectors (i.e. norm 0) are kept empty.
af::array normalize_handle_zero_vectors(const af::array& a) {
    af::array norm_a = af::tile(vecnorm(a), 1, 1, 1, 3);
    af::array normalized = a / norm_a;
    af::replace(normalized, norm_a != 0, 0);
    return normalized;
}

af::array vecnorm(const af::array& a) { return af::sqrt(af::sum(a * a, 3)); }

// Mean value of i = 0, 1, 2 entry entry
double meani(const af::array& a, const int i) {
    return af::mean(af::mean(af::mean(a(af::span, af::span, af::span, i), 0), 1), 2).scalar<double>();
}

// Frobenius Norm
//||A||=sqrt(sum(fabs(a)))
double FrobeniusNorm(const af::array& a) {
    return af::sqrt(af::mean(af::mean(af::mean(af::mean(a * a, 0), 1), 2), 3)).scalar<double>();
}

// Experimental: eucledian norm
double euclnorm(const af::array& a) {
    return af::mean(af::mean(af::mean(af::mean((a * a), 0), 1), 2), 3).scalar<double>();
}

/// Mean of absolute difference
double mean_abs_diff(const af::array& a, const af::array& b) {
    return afvalue_as_f64(af::mean(af::mean(af::mean(af::mean(af::abs(a - b), 0), 1), 2), 3));
}

/// Mean of relative difference
double mean_rel_diff(const af::array& first, const af::array& second) {
    af::array temp = af::abs(2 * (first - second) / (first + second));
    af::replace(temp, first != 0 || second != 0,
                af::constant(0., temp.dims(),
                             f64)); // Avoiding division by zero: setting element to
                                    // zero if both input elements are zero
    return afvalue_as_f64(af::mean(af::mean(af::mean(af::mean(temp, 0), 1), 2), 3));
}

/// Max of absolute difference
double max_abs_diff(const af::array& a, const af::array& b) {
    return afvalue_as_f64(af::max(af::max(af::max(af::max(af::abs(a - b), 0), 1), 2), 3));
}

/// Max of relative difference
double max_rel_diff(const af::array& first, const af::array& second) {
    af::array temp = af::abs(2 * (first - second) / (first + second));
    af::replace(temp, first != 0 || second != 0,
                af::constant(0., temp.dims(),
                             f64)); // Avoiding division by zero: setting element to
                                    // zero if both input elements are zero
    return afvalue_as_f64(af::max(af::max(af::max(af::max(temp, 0), 1), 2), 3));
}

/// Absolute difference less than precision: Element-wise comparision of
/// absolute difference of two arrays. Checks whether | x - y | < precision.
/// Returns true if all values are below precision and false otherwise.
bool abs_diff_lt_precision(const af::array& first, const af::array& second, double precision, bool verbose) {
    unsigned int zero_if_equal =
        afvalue_u32(af::sum(af::sum(af::sum(af::sum(!(af::abs(first - second) < precision), 0), 1), 2), 3));
    if (verbose) {
        if (zero_if_equal == 0) {
            std::cout << "\33[1;32mSucess:\33[0m All " << first.elements()
                      << " absolute values of element-wise differences are "
                         "below precision of "
                      << precision << std::endl;
        } else {
            std::cout << "\33[1;31mError!\33[0m " << zero_if_equal << " out of " << first.elements()
                      << " absolute values of element-wise differences are "
                         "above precision of "
                      << precision << std::endl;
        }
    }
    if (zero_if_equal == 0) {
        return true;
    } else {
        return false;
    }
}

/// Relative difference less than precision: Element-wise comparision of
/// relative difference of two arrays. Checks whether | 2(x-y)/(x+y) | <
/// precision. Returns true if all values are below precision and false
/// otherwise.
bool rel_diff_lt_precision(const af::array& first, const af::array& second, double precision, bool verbose) {
    af::array temp = af::abs(2 * (first - second) / (first + second));
    af::replace(temp, first != 0 || second != 0,
                af::constant(0., temp.dims(),
                             f64)); // set element to zero if both input elements are zero
    unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum(!(temp < precision), 0), 1), 2), 3));
    if (verbose) {
        if (zero_if_equal == 0) {
            std::cout << "\33[1;32mSucess:\33[0m All " << first.elements()
                      << " relative values of element-wise differences are "
                         "below precision of "
                      << precision << std::endl;
        } else {
            std::cout << "\33[1;31mError!\33[0m " << zero_if_equal << " out of " << first.elements()
                      << " relative values of element-wise differences are "
                         "above precision of "
                      << precision << std::endl;
        }
    }
    if (zero_if_equal == 0) {
        return true;
    } else {
        return false;
    }
}

/// Upper bound for absolute difference
double abs_diff_upperbound(const af::array& a, const af::array& b, bool verbose, double start_precision, double factor1,
                           double factor2) {
    double prec = start_precision;
    double prec_prev = prec;
    while (abs_diff_lt_precision(a, b, prec, false) and prec > 1e-300) {
        if (verbose) {
            std::cout << "prec = " << prec << std::endl;
        }
        prec_prev = prec;
        prec = factor1 * prec;
    }

    prec = prec_prev;
    while (abs_diff_lt_precision(a, b, prec, false) and prec > 1e-300) {
        if (verbose) {
            std::cout << "prec = " << prec << std::endl;
        }
        prec_prev = prec;
        prec = factor2 * prec;
    }
    return prec_prev;
}

/// Upper bound for relative difference
double rel_diff_upperbound(const af::array& a, const af::array& b, bool verbose, double start_precision, double factor1,
                           double factor2) {
    double prec = start_precision;
    double prec_prev = prec;
    while (rel_diff_lt_precision(a, b, prec, false) and prec > 1e-300) {
        if (verbose) {
            std::cout << "prec = " << prec << std::endl;
        }
        prec_prev = prec;
        prec = factor1 * prec;
    }

    prec = prec_prev;
    while (rel_diff_lt_precision(a, b, prec, false) and prec > 1e-300) {
        if (verbose) {
            std::cout << "prec = " << prec << std::endl;
        }
        prec_prev = prec;
        prec = factor2 * prec;
    }
    return prec_prev;
}

// Minimum value
double min_4d(const af::array& a) { return af::min(af::min(af::min(af::min(a.as(f64), 0), 1), 2), 3).scalar<double>(); }

std::pair<int, int> k2ij(const int k, const int n) {
    const int i = n - 1 - static_cast<int>(std::floor(std::sqrt(-8 * k + 4 * n * (n + 1) - 7) / 2.0 - 0.5));
    const int j = k + i - n * (n + 1) / 2 + (n - i) * ((n - i) + 1) / 2;
    return std::make_pair(i, j);
}

int ij2k(const int i, const int j, const int n) { return (n * (n + 1) / 2) - (n - i) * ((n - i) + 1) / 2 + j - i; }

unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int ni) { return i + ni * j; }
unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int ni,
                    const unsigned int nj) {
    return i + ni * (j + nj * k);
}
unsigned int stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l,
                    const unsigned int ni, const unsigned int nj, const unsigned int nk) {
    return i + ni * (j + nj * (k + nk * l));
}

af::randomEngine rand_engine_current_time() {
    unsigned long long int seed = static_cast<unsigned long long int>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
            .count());
    return af::randomEngine(AF_RANDOM_ENGINE_DEFAULT, seed);
}

// TODO check with c++14 (we used uncommented due to incompability with c++11
// needed by cython)
////RK4 based on https://rosettacode.org/wiki/Runge-Kutta_method
// auto rk4(af::array f(double, af::array))
//{
//        return
//        [       f            ](double t,  af::array y, double dt ) ->
//        af::array { return [t, y, dt, f            ]( af::array  dy1) ->
//        af::array { return [t, y, dt, f, dy1        ]( af::array  dy2) ->
//        af::array { return [t, y, dt, f, dy1, dy2    ]( af::array  dy3) ->
//        af::array { return [t, y, dt, f, dy1, dy2, dy3]( af::array  dy4) ->
//        af::array { return ( dy1 + 2*dy2 + 2*dy3 + dy4 ) / 6   ;} ( dt * f(
//        t+dt  , y+dy3   )          );} ( dt * f( t+dt/2, y+dy2/2 ) );} ( dt *
//        f( t+dt/2, y+dy1/2 )          );} ( dt * f( t     , y       ) );} ;
//}
// TODO END

// int main(void)
//{
//        const double TIME_MAXIMUM = 10.0, WHOLE_TOLERANCE = 1e-12 ;
//        const double T_START = 0.0, Y_START = 1.0, DT = 0.10;
//
//        auto eval_diff_eqn = [               ](double t, double y)->double{
//        return t*sqrt(y)                         ; } ; auto eval_solution = [
//        ](double t          )->double{ return pow(t*t+4, 2)/16 ; } ; auto
//        find_error    = [eval_solution  ](double t, double y)->double{ return
//        fabs(y-eval_solution(t))          ; } ; auto is_whole      =
//        [WHOLE_TOLERANCE](double t          )->bool  { return fabs(t-round(t))
//        < WHOLE_TOLERANCE; } ;
//
//        auto dy = rk4( eval_diff_eqn ) ;
//
//        double y = Y_START, t = T_START ;
//
//        while(t <= TIME_MAXIMUM) {
//          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n", t,
//          y, find_error(t, y)); } y += dy(t, y, DT) ; t += DT;
//        }
//        return 0;
//}

} // namespace magnumaf::util
