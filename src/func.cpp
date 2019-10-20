#include "func.hpp"

namespace magnumaf{

WrappedArray::WrappedArray(af::array array) : array(array) {
}

WrappedArray::WrappedArray(long int  array_ptr){
    set_array(array_ptr);
}

void WrappedArray::set_array(long int array_ptr){
    void **a = (void **)array_ptr;
    this->array = *( new af::array( *a ));
}

long int WrappedArray::get_array_addr(){
    af::array *a = new af::array(this->array);
    return (long int) a->get();
}

af::array cross4(const af::array& a, const af::array& b){
  af::array c= af::array(a.dims(0), a.dims(1), a.dims(2), 3, f32);
  c(af::span, af::span, af::span, 0)=a(af::span, af::span, af::span, 1)*b(af::span, af::span, af::span, 2)-a(af::span, af::span, af::span, 2)*b(af::span, af::span, af::span, 1);
  c(af::span, af::span, af::span, 1)=a(af::span, af::span, af::span, 2)*b(af::span, af::span, af::span, 0)-a(af::span, af::span, af::span, 0)*b(af::span, af::span, af::span, 2);
  c(af::span, af::span, af::span, 2)=a(af::span, af::span, af::span, 0)*b(af::span, af::span, af::span, 1)-a(af::span, af::span, af::span, 1)*b(af::span, af::span, af::span, 0);
  return c;
};

af::array dotproduct(const af::array& a, const af::array& b){
  return sum(a*b, 3);
}

/// Returns the value of array with only one element
float afvalue(const af::array& a){
    if (a.dims(0) != 1 || a.dims(1) != 1 || a.dims(2) != 1 || a.dims(3) != 1){
        std::cout << "\33[1;31mWarning:\33[0m afvalue requested from array with dim4 =[" << a.dims() <<"] != [1 1 1 1]. Only first entry will be returned. This may lead to unexpected behaviour." << std::endl;
    }
    float *dhost=NULL;
    dhost = a.host<float>();
    float value = dhost[0];
    af::freeHost(dhost);
    return value;
}

unsigned int afvalue_u32(const af::array& a){
  unsigned int *dhost=NULL;
  dhost = a.host<unsigned int>();
  unsigned int value = dhost[0];
  af::freeHost(dhost);
  return value;
}

float full_inner_product(const af::array& a, const af::array& b){
  return afvalue(sum(sum(sum(sum(a*b, 3), 2), 1), 0));
}

af::array renormalize(const af::array& a){
  return a/tile(sqrt(sum(a*a, 3)), 1, 1, 1, 3);
}

//Renormalization where values with Ms zero are set from inf to zero
af::array renormalize_handle_zero_values(const af::array& a){
    af::array norm_a = tile(sqrt(sum(a*a, 3)), 1, 1, 1, 3);
    af::array normalized = a/norm_a;
    replace(normalized, norm_a != 0, 0);
    return normalized;

    //TODO for (af::array& a) only: return replace(a/tile(sqrt(sum(a*a, 3)), 1, 1, 1, 3), a!=0., 0.);
}

af::array vecnorm(const af::array& a){
  return sqrt(sum(a*a, 3));
}


//Mean value of i = 0, 1, 2 entry entry
float meani(const af::array& a, const int i){
  float *norm_host=NULL;
  norm_host = mean(mean(mean(a(af::span, af::span, af::span, i), 0), 1), 2).host<float>();
  float norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}
//Frobenius Norm
//||A||=sqrt(sum(fabs(a)))
float FrobeniusNorm(const af::array& a){
  float *norm_host=NULL;
  norm_host = sqrt(mean(mean(mean(mean(a*a, 0), 1), 2), 3)).host<float>();
  float norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}

//Experimental: eucledian norm
float euclnorm(const af::array& a){
  float *norm_host=NULL;
  norm_host = mean(mean(mean(mean((a*a), 0), 1), 2), 3).host<float>();
  float norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}

/// Mean of absolute difference
float mean_abs_diff(const af::array& a, const af::array& b){
    return afvalue(af::mean(af::mean(af::mean(af::mean(af::abs(a - b), 0), 1), 2), 3));
}

/// Mean of relative difference
float mean_rel_diff(const af::array& first, const af::array& second){
    af::array temp = af::abs(2*(first - second)/(first + second));
    af::replace(temp, first!=0 || second!=0, af::constant(0., temp.dims(), f32));// Avoiding division by zero: setting element to zero if both input elements are zero
    return afvalue(af::mean(af::mean(af::mean(af::mean(temp, 0), 1), 2), 3));
}

/// Max of absolute difference
float max_abs_diff(const af::array& a, const af::array& b){
    return afvalue(af::max(af::max(af::max(af::max(af::abs(a - b), 0), 1), 2), 3));
}

/// Max of relative difference
float max_rel_diff(const af::array& first, const af::array& second){
    af::array temp = af::abs(2*(first - second)/(first + second));
    af::replace(temp, first!=0 || second!=0, af::constant(0., temp.dims(), f32));// Avoiding division by zero: setting element to zero if both input elements are zero
    return afvalue(af::max(af::max(af::max(af::max(temp, 0), 1), 2), 3));
}

/// Absolute difference less than precision: Element-wise comparision of absolute difference of two arrays. Checks whether | x - y | < precision. Returns true if all values are below precision and false otherwise.
bool abs_diff_lt_precision(af::array first, af::array second, float precision, bool verbose){
    unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum( !(af::abs(first - second) < precision), 0), 1), 2), 3));
    if (verbose){
        if (zero_if_equal == 0) std::cout << "\33[1;32mSucess:\33[0m All " << first.elements() << " absolute values of element-wise differences are below precision of " << precision << std::endl;
        else {
            std::cout << "\33[1;31mError!\33[0m " << zero_if_equal << " out of " << first.elements() << " absolute values of element-wise differences are above precision of " << precision << std::endl;
        }
    }
    if (zero_if_equal == 0) return true;
    else return false;
}

/// Relative difference less than precision: Element-wise comparision of relative difference of two arrays. Checks whether | 2(x-y)/(x+y) | < precision. Returns true if all values are below precision and false otherwise.
bool rel_diff_lt_precision(af::array first, af::array second, float precision, bool verbose){
    af::array temp = af::abs(2*(first - second)/(first + second));
    af::replace(temp, first!=0 || second!=0, af::constant(0., temp.dims(), f32));//set element to zero if both input elements are zero
    unsigned int zero_if_equal = afvalue_u32(af::sum(af::sum(af::sum(af::sum( !(temp < precision), 0), 1), 2), 3));
    if (verbose){
        if (zero_if_equal == 0) std::cout << "\33[1;32mSucess:\33[0m All " << first.elements() << " relative values of element-wise differences are below precision of " << precision << std::endl;
        else {
            std::cout << "\33[1;31mError!\33[0m " << zero_if_equal << " out of " << first.elements() << " relative values of element-wise differences are above precision of " << precision << std::endl;
        }
    }
    if (zero_if_equal == 0) return true;
    else return false;
}

/// Upper bound for absolute difference
float abs_diff_upperbound(const af::array& a, const af::array& b, bool verbose, float start_precision, float factor1, float factor2){
    float prec = start_precision;
    float prec_prev = prec;
    while(abs_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose) std::cout << "prec = " <<prec << std::endl;
        prec_prev = prec;
        prec = factor1 * prec;
    }

    prec = prec_prev;
    while(abs_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose) std::cout << "prec = " <<prec << std::endl;
        prec_prev = prec;
        prec = factor2 * prec;
    }
    return prec_prev;
}

/// Upper bound for relative difference
float rel_diff_upperbound(const af::array& a, const af::array& b, bool verbose, float start_precision, float factor1, float factor2){
    float prec = start_precision;
    float prec_prev = prec;
    while(rel_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose) std::cout << "prec = " <<prec << std::endl;
        prec_prev = prec;
        prec = factor1 * prec;
    }

    prec = prec_prev;
    while(rel_diff_lt_precision(a, b, prec, false) and prec > 1e-300)
    {
        if (verbose) std::cout << "prec = " <<prec << std::endl;
        prec_prev = prec;
        prec = factor2 * prec;
    }
    return prec_prev;
}

//Experimental: eucledian norm
//float maxnorm(const af::array& a){
//  float *maxnorm_host=NULL;
//  maxnorm_host = mean(mean(mean(mean((a*a), 0), 1), 2), 3).host<float>();
//  //maxnorm_host = max(max(max(max(abs(a), 0), 1), 2), 3).host<float>();
//  float maxnorm = maxnorm_host[0];
//  af::freeHost(maxnorm_host);
//  return maxnorm;
//}

//void calcm(State state, LLG Llg, std::ostream& myfile){
//  float *host_mx=NULL, *host_my=NULL, *host_mz=NULL;
//  host_mx = mean(mean(mean(state.m(af::span, af::span, af::span, 0), 0), 1), 2).host<float>();
//  host_my = mean(mean(mean(state.m(af::span, af::span, af::span, 1), 0), 1), 2).host<float>();
//  host_mz = mean(mean(mean(state.m(af::span, af::span, af::span, 2), 0), 1), 2).host<float>();
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
float maxnorm(const af::array& a){
  float *maxnorm_host=NULL;
  maxnorm_host = max(max(max(max(abs(a), 0), 1), 2), 3).host<float>();
  float maxnorm = maxnorm_host[0];
  af::freeHost(maxnorm_host);
  return maxnorm;
}

// Minimum value
float minval(const af::array& a){
  float *minval_host=NULL;
  minval_host = min(min(min(min(a, 0), 1), 2), 3).host<float>();
  float minval = minval_host[0];
  af::freeHost(minval_host);
  return minval;
}

std::pair<int, int> util::k2ij(const int k, const int n){
    const int i = n - 1 - std::floor(std::sqrt(-8*k + 4*n*(n+1)-7)/2.0 - 0.5);
    const int j = k + i - n*(n+1)/2 + (n-i)*((n-i)+1)/2;
    return std::make_pair(i, j);
}

int util::ij2k(const int i, const int j, const int n){
    return (n*(n+1)/2) - (n-i)*((n-i)+1)/2 + j - i;
}

unsigned int util::stride(const unsigned int i, const unsigned int j, const unsigned int ni){return i + ni * j;}
unsigned int util::stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int ni, const unsigned int nj){return i + ni * (j + nj * k);}
unsigned int util::stride(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l, const unsigned int ni, const unsigned int nj, const unsigned int nk){return i + ni * (j + nj * (k + nk * l));}



af::randomEngine util::rand_engine_current_time(){
    unsigned long long int seed = std::chrono::duration_cast< std::chrono::nanoseconds >\
                                  (std::chrono::system_clock::now().time_since_epoch()).count();
    return af::randomEngine(AF_RANDOM_ENGINE_DEFAULT, seed);
}



//TODO check with c++14 (we used uncommented due to incompability with c++11 needed by cython)
////RK4 based on https://rosettacode.org/wiki/Runge-Kutta_method
//auto rk4(af::array f(float, af::array))
//{
//        return
//        [       f            ](float t,  af::array y, float dt ) -> af::array { return
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
//        const float TIME_MAXIMUM = 10.0, WHOLE_TOLERANCE = 1e-12 ;
//        const float T_START = 0.0, Y_START = 1.0, DT = 0.10;
//
//        auto eval_diff_eqn = [               ](float t, float y)->float{ return t*sqrt(y)                         ; } ;
//        auto eval_solution = [               ](float t          )->float{ return pow(t*t+4, 2)/16                   ; } ;
//        auto find_error    = [eval_solution  ](float t, float y)->float{ return fabs(y-eval_solution(t))          ; } ;
//        auto is_whole      = [WHOLE_TOLERANCE](float t          )->bool  { return fabs(t-round(t)) < WHOLE_TOLERANCE; } ;
//
//        auto dy = rk4( eval_diff_eqn ) ;
//
//        float y = Y_START, t = T_START ;
//
//        while(t <= TIME_MAXIMUM) {
//          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n", t, y, find_error(t, y)); }
//          y += dy(t, y, DT) ; t += DT;
//        }
//        return 0;
//}

}// namespace magnumaf
