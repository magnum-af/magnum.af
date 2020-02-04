#include "func.hpp"

namespace magnumafcpp
{


std::array<double, 3> spacial_mean_in_region(const af::array& vectorfield, const af::array& region){
    if( vectorfield.dims(3) != 3 or region.dims(3) != 1)
    {
        printf("\33[1;31mWarning:\33[0m spacial_mean_in_region: dimensions do not match, returning [-1, -1, -1]\n");
        return {-1, -1, -1};
    }
    else{
        af::array zero_if_input_is_zero_else_one = ! af::iszero(region);
        const unsigned number_of_nonzero_elements = af::sum(af::sum(af::sum(zero_if_input_is_zero_else_one, 0), 1), 2).scalar<unsigned>();
        af::array considered_values = vectorfield * af::tile(zero_if_input_is_zero_else_one, 1, 1, 1, 3);
        af::array sum_considered_values = af::sum(af::sum(af::sum(considered_values, 0), 1), 2);
        af::array mean_considered_values = sum_considered_values/number_of_nonzero_elements;
        double xmean = mean_considered_values(0, 0, 0, 0).scalar<double>();
        double ymean = mean_considered_values(0, 0, 0, 1).scalar<double>();
        double zmean = mean_considered_values(0, 0, 0, 2).scalar<double>();
        return {xmean, ymean, zmean};
    }
}

WrappedArray::WrappedArray(af::array array) : array(array)
{
}

WrappedArray::WrappedArray(long int array_ptr)
{
    set_array(array_ptr);
}

void WrappedArray::set_array(long int array_ptr)
{
    void **a = (void **)array_ptr;
    this->array = *(new af::array(*a));
}

long int WrappedArray::get_array_addr()
{
    af::array *a = new af::array(this->array);
    return (long int)a->get();
}

af::array cross4(const af::array &a, const af::array &b)
{
    af::array c = af::array(a.dims(0), a.dims(1), a.dims(2), 3, f64);
    c(af::span, af::span, af::span, 0) = a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 2) - a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 1);
    c(af::span, af::span, af::span, 1) = a(af::span, af::span, af::span, 2) * b(af::span, af::span, af::span, 0) - a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 2);
    c(af::span, af::span, af::span, 2) = a(af::span, af::span, af::span, 0) * b(af::span, af::span, af::span, 1) - a(af::span, af::span, af::span, 1) * b(af::span, af::span, af::span, 0);
    return c;
}

/// Cross product for vector fields of format [nx, ny, nz, 3]
af::array cross4shift(const af::array &a, const af::array &b)
{
    af::array ashift = af::shift(a, 0, 0, 0, -1);
    af::array ashift2 = af::shift(a, 0, 0, 0, -2);
    af::array bshift = af::shift(b, 0, 0, 0, -2);
    af::array bshift2 = af::shift(b, 0, 0, 0, -1);
    return ashift * bshift - ashift2 * bshift2;
}

af::array dotproduct(const af::array &a, const af::array &b)
{
    return sum(a * b, 3);
}

/// Returns the value of array with only one element
double afvalue(const af::array &a)
{
    if (a.dims(0) != 1 || a.dims(1) != 1 || a.dims(2) != 1 || a.dims(3) != 1)
    {
        std::cout << "\33[1;31mWarning:\33[0m afvalue requested from array with dim4 =[" << a.dims() << "] != [1 1 1 1]. Only first entry will be returned. This may lead to unexpected behaviour." << std::endl;
    }
    double *dhost = NULL;
    dhost = a.host<double>();
    double value = dhost[0];
    af::freeHost(dhost);
    return value;
}

unsigned int afvalue_u32(const af::array &a)
{
    unsigned int *dhost = NULL;
    dhost = a.host<unsigned int>();
    unsigned int value = dhost[0];
    af::freeHost(dhost);
    return value;
}

double full_inner_product(const af::array &a, const af::array &b)
{
    return afvalue(sum(sum(sum(sum(a * b, 3), 2), 1), 0));
}

af::array renormalize(const af::array &a)
{
    return a / tile(sqrt(sum(a * a, 3)), 1, 1, 1, 3);
}

//Renormalization where values with Ms zero are set from inf to zero
af::array renormalize_handle_zero_values(const af::array &a)
{
    af::array norm_a = tile(sqrt(sum(a * a, 3)), 1, 1, 1, 3);
    af::array normalized = a / norm_a;
    replace(normalized, norm_a != 0, 0);
    return normalized;

    //TODO for (af::array& a) only: return replace(a/tile(sqrt(sum(a*a, 3)), 1, 1, 1, 3), a!=0., 0.);
}

af::array vecnorm(const af::array &a)
{
    return sqrt(sum(a * a, 3));
}

//Mean value of i = 0, 1, 2 entry entry
double meani(const af::array &a, const int i)
{
    double *norm_host = NULL;
    norm_host = mean(mean(mean(a(af::span, af::span, af::span, i), 0), 1), 2).host<double>();
    double norm = norm_host[0];
    af::freeHost(norm_host);
    return norm;
}
//Frobenius Norm
//||A||=sqrt(sum(fabs(a)))
double FrobeniusNorm(const af::array &a)
{
    double *norm_host = NULL;
    norm_host = sqrt(mean(mean(mean(mean(a * a, 0), 1), 2), 3)).host<double>();
    double norm = norm_host[0];
    af::freeHost(norm_host);
    return norm;
}

//Experimental: eucledian norm
double euclnorm(const af::array &a)
{
    double *norm_host = NULL;
    norm_host = mean(mean(mean(mean((a * a), 0), 1), 2), 3).host<double>();
    double norm = norm_host[0];
    af::freeHost(norm_host);
    return norm;
}

/// Mean of absolute difference
double mean_abs_diff(const af::array &a, const af::array &b)
{
    return afvalue(af::mean(af::mean(af::mean(af::mean(af::abs(a - b), 0), 1), 2), 3));
}

/// Mean of relative difference
double mean_rel_diff(const af::array &first, const af::array &second)
{
    af::array temp = af::abs(2 * (first - second) / (first + second));
    af::replace(temp, first != 0 || second != 0, af::constant(0., temp.dims(), f64)); // Avoiding division by zero: setting element to zero if both input elements are zero
    return afvalue(af::mean(af::mean(af::mean(af::mean(temp, 0), 1), 2), 3));
}

/// Max of absolute difference
double max_abs_diff(const af::array &a, const af::array &b)
{
    return afvalue(af::max(af::max(af::max(af::max(af::abs(a - b), 0), 1), 2), 3));
}

/// Max of relative difference
double max_rel_diff(const af::array &first, const af::array &second)
{
    af::array temp = af::abs(2 * (first - second) / (first + second));
    af::replace(temp, first != 0 || second != 0, af::constant(0., temp.dims(), f64)); // Avoiding division by zero: setting element to zero if both input elements are zero
    return afvalue(af::max(af::max(af::max(af::max(temp, 0), 1), 2), 3));
}

/// Absolute difference less than precision: Element-wise comparision of absolute difference of two arrays. Checks whether | x - y | < precision. Returns true if all values are below precision and false otherwise.
bool abs_diff_lt_precision(af::array first, af::array second, double precision, bool verbose)
{
    unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum(!(af::abs(first - second) < precision), 0), 1), 2), 3));
    if (verbose)
    {
        if (zero_if_equal == 0)
            std::cout << "\33[1;32mSucess:\33[0m All " << first.elements() << " absolute values of element-wise differences are below precision of " << precision << std::endl;
        else
        {
            std::cout << "\33[1;31mError!\33[0m " << zero_if_equal << " out of " << first.elements() << " absolute values of element-wise differences are above precision of " << precision << std::endl;
        }
    }
    if (zero_if_equal == 0)
        return true;
    else
        return false;
}

/// Relative difference less than precision: Element-wise comparision of relative difference of two arrays. Checks whether | 2(x-y)/(x+y) | < precision. Returns true if all values are below precision and false otherwise.
bool rel_diff_lt_precision(af::array first, af::array second, double precision, bool verbose)
{
    af::array temp = af::abs(2 * (first - second) / (first + second));
    af::replace(temp, first != 0 || second != 0, af::constant(0., temp.dims(), f64)); //set element to zero if both input elements are zero
    unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum(!(temp < precision), 0), 1), 2), 3));
    if (verbose)
    {
        if (zero_if_equal == 0)
            std::cout << "\33[1;32mSucess:\33[0m All " << first.elements() << " relative values of element-wise differences are below precision of " << precision << std::endl;
        else
        {
            std::cout << "\33[1;31mError!\33[0m " << zero_if_equal << " out of " << first.elements() << " relative values of element-wise differences are above precision of " << precision << std::endl;
        }
    }
    if (zero_if_equal == 0)
        return true;
    else
        return false;
}

/// Upper bound for absolute difference
double abs_diff_upperbound(const af::array &a, const af::array &b, bool verbose, double start_precision, double factor1, double factor2)
{
    double prec = start_precision;
    double prec_prev = prec;
    while (abs_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose)
            std::cout << "prec = " << prec << std::endl;
        prec_prev = prec;
        prec = factor1 * prec;
    }

    prec = prec_prev;
    while (abs_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose)
            std::cout << "prec = " << prec << std::endl;
        prec_prev = prec;
        prec = factor2 * prec;
    }
    return prec_prev;
}

/// Upper bound for relative difference
double rel_diff_upperbound(const af::array &a, const af::array &b, bool verbose, double start_precision, double factor1, double factor2)
{
    double prec = start_precision;
    double prec_prev = prec;
    while (rel_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose)
            std::cout << "prec = " << prec << std::endl;
        prec_prev = prec;
        prec = factor1 * prec;
    }

    prec = prec_prev;
    while (rel_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose)
            std::cout << "prec = " << prec << std::endl;
        prec_prev = prec;
        prec = factor2 * prec;
    }
    return prec_prev;
}

//Experimental: eucledian norm
//double maxnorm(const af::array& a){
//  double *maxnorm_host=NULL;
//  maxnorm_host = mean(mean(mean(mean((a*a), 0), 1), 2), 3).host<double>();
//  //maxnorm_host = max(max(max(max(abs(a), 0), 1), 2), 3).host<double>();
//  double maxnorm = maxnorm_host[0];
//  af::freeHost(maxnorm_host);
//  return maxnorm;
//}

//void calcm(State state, LLG Llg, std::ostream& myfile){
//  double *host_mx=NULL, *host_my=NULL, *host_mz=NULL;
//  host_mx = mean(mean(mean(state.m(af::span, af::span, af::span, 0), 0), 1), 2).host<double>();
//  host_my = mean(mean(mean(state.m(af::span, af::span, af::span, 1), 0), 1), 2).host<double>();
//  host_mz = mean(mean(mean(state.m(af::span, af::span, af::span, 2), 0), 1), 2).host<double>();
//  myfile << std::setw(12) << state.t << "\t";
//  myfile << std::setw(12) << state.t*1e9 << "\t" << host_mx[0] << "\t"<< host_my[0] << "\t";
//  myfile << host_mz[0] <<"\t" << Llg.E(state) <<std::endl;
//
//  af::freeHost(host_mx);
//  af::freeHost(host_my);
//  af::freeHost(host_mz);
//
//}

// Maximum value norm
double maxnorm(const af::array &a)
{
    double *maxnorm_host = NULL;
    maxnorm_host = max(max(max(max(abs(a), 0), 1), 2), 3).host<double>();
    double maxnorm = maxnorm_host[0];
    af::freeHost(maxnorm_host);
    return maxnorm;
}

// Minimum value
double minval(const af::array &a)
{
    double *minval_host = NULL;
    minval_host = min(min(min(min(a, 0), 1), 2), 3).host<double>();
    double minval = minval_host[0];
    af::freeHost(minval_host);
    return minval;
}

std::pair<int, int> util::k2ij(const int k, const int n)
{
    const int i = n - 1 - static_cast<int>(std::floor(std::sqrt(-8 * k + 4 * n * (n + 1) - 7) / 2.0 - 0.5));
    const int j = k + i - n * (n + 1) / 2 + (n - i) * ((n - i) + 1) / 2;
    return std::make_pair(i, j);
}

int util::ij2k(const int i, const int j, const int n)
{
    return (n * (n + 1) / 2) - (n - i) * ((n - i) + 1) / 2 + j - i;
}

unsigned int util::stride(const unsigned int i, const unsigned int j, const unsigned int ni) { return i + ni * j; }
unsigned int util::stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int ni, const unsigned int nj) { return i + ni * (j + nj * k); }
unsigned int util::stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l, const unsigned int ni, const unsigned int nj, const unsigned int nk) { return i + ni * (j + nj * (k + nk * l)); }

af::randomEngine util::rand_engine_current_time()
{
    unsigned long long int seed = static_cast<unsigned long long int>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    return af::randomEngine(AF_RANDOM_ENGINE_DEFAULT, seed);
}

//TODO check with c++14 (we used uncommented due to incompability with c++11 needed by cython)
////RK4 based on https://rosettacode.org/wiki/Runge-Kutta_method
//auto rk4(af::array f(double, af::array))
//{
//        return
//        [       f            ](double t,  af::array y, double dt ) -> af::array { return
//        [t, y, dt, f            ](                    af::array  dy1) -> af::array { return
//        [t, y, dt, f, dy1        ](                    af::array  dy2) -> af::array { return
//        [t, y, dt, f, dy1, dy2    ](                    af::array  dy3) -> af::array { return
//        [t, y, dt, f, dy1, dy2, dy3](                    af::array  dy4) -> af::array { return
//        ( dy1 + 2*dy2 + 2*dy3 + dy4 ) / 6   ;} (
//        dt * f( t+dt  , y+dy3   )          );} (
//        dt * f( t+dt/2, y+dy2/2 )          );} (
//        dt * f( t+dt/2, y+dy1/2 )          );} (
//        dt * f( t     , y       )          );} ;
//}
//TODO END

//int main(void)
//{
//        const double TIME_MAXIMUM = 10.0, WHOLE_TOLERANCE = 1e-12 ;
//        const double T_START = 0.0, Y_START = 1.0, DT = 0.10;
//
//        auto eval_diff_eqn = [               ](double t, double y)->double{ return t*sqrt(y)                         ; } ;
//        auto eval_solution = [               ](double t          )->double{ return pow(t*t+4, 2)/16                   ; } ;
//        auto find_error    = [eval_solution  ](double t, double y)->double{ return fabs(y-eval_solution(t))          ; } ;
//        auto is_whole      = [WHOLE_TOLERANCE](double t          )->bool  { return fabs(t-round(t)) < WHOLE_TOLERANCE; } ;
//
//        auto dy = rk4( eval_diff_eqn ) ;
//
//        double y = Y_START, t = T_START ;
//
//        while(t <= TIME_MAXIMUM) {
//          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n", t, y, find_error(t, y)); }
//          y += dy(t, y, DT) ; t += DT;
//        }
//        return 0;
//}

} // namespace magnumafcpp
