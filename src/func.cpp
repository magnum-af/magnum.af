#include "func.hpp"
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

af::array cross4(const af::array& a,const af::array& b){
  af::array c= af::array(a.dims(0),a.dims(1),a.dims(2),3,f64);
  c(af::span,af::span,af::span,0)=a(af::span,af::span,af::span,1)*b(af::span,af::span,af::span,2)-a(af::span,af::span,af::span,2)*b(af::span,af::span,af::span,1);
  c(af::span,af::span,af::span,1)=a(af::span,af::span,af::span,2)*b(af::span,af::span,af::span,0)-a(af::span,af::span,af::span,0)*b(af::span,af::span,af::span,2);
  c(af::span,af::span,af::span,2)=a(af::span,af::span,af::span,0)*b(af::span,af::span,af::span,1)-a(af::span,af::span,af::span,1)*b(af::span,af::span,af::span,0);
  return c;
};

af::array dotproduct(const af::array& a, const af::array& b){
  return sum(a*b,3);
}

double afvalue(const af::array& a){
  double *dhost=NULL;
  dhost = a.host<double>();
  double value = dhost[0];
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

double full_inner_product(const af::array& a, const af::array& b){
  return afvalue(sum(sum(sum(sum(a*b,3),2),1),0));
}

af::array renormalize(const af::array& a){
  return a/tile(sqrt(sum(a*a,3)),1,1,1,3); 
}

//Renormalization where values with Ms zero are set from inf to zero
af::array renormalize_handle_zero_values(const af::array& a){
    af::array norm_a = tile(sqrt(sum(a*a,3)),1,1,1,3); 
    af::array renorm = a/norm_a;
    replace(renorm,norm_a!=0,0);
    return renorm;
    
    //TODO for (af::array& a) only: return replace(a/tile(sqrt(sum(a*a,3)),1,1,1,3), a!=0., 0.);
}

af::array vecnorm(const af::array& a){
  return sqrt(sum(a*a,3)); 
}


//Mean value of i = 0,1,2 entry entry
double meani(const af::array& a, const int i){
  double *norm_host=NULL;
  norm_host = mean(mean(mean(a(af::span,af::span,af::span,i),0),1),2).host<double>();
  double norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}
//Frobenius Norm
//||A||=sqrt(sum(fabs(a)))
double FrobeniusNorm(const af::array& a){
  double *norm_host=NULL;
  norm_host = sqrt(mean(mean(mean(mean(a*a,0),1),2),3)).host<double>();
  double norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}

//Experimental: eucledian norm
double euclnorm(const af::array& a){
  double *norm_host=NULL;
  norm_host = mean(mean(mean(mean((a*a),0),1),2),3).host<double>();
  double norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}


/// Absolute difference less than precision: Element-wise comparision of absolute difference of two arrays. Checks whether | x - y | < precision. Returns true if all values are below precision and false otherwise.
bool abs_diff_lt_precision(af::array first, af::array second, double precision, bool verbose){
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

//Experimental: eucledian norm
//double maxnorm(const af::array& a){
//  double *maxnorm_host=NULL;
//  maxnorm_host = mean(mean(mean(mean((a*a),0),1),2),3).host<double>();
//  //maxnorm_host = max(max(max(max(abs(a),0),1),2),3).host<double>();
//  double maxnorm = maxnorm_host[0];
//  af::freeHost(maxnorm_host);
//  return maxnorm;
//}

//void calcm(State state, LLG Llg, std::ostream& myfile){
//  double *host_mx=NULL, *host_my=NULL, *host_mz=NULL;
//  host_mx = mean(mean(mean(state.m(af::span,af::span,af::span,0),0),1),2).host<double>();
//  host_my = mean(mean(mean(state.m(af::span,af::span,af::span,1),0),1),2).host<double>();
//  host_mz = mean(mean(mean(state.m(af::span,af::span,af::span,2),0),1),2).host<double>();
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
double maxnorm(const af::array& a){
  double *maxnorm_host=NULL;
  maxnorm_host = max(max(max(max(abs(a),0),1),2),3).host<double>();
  double maxnorm = maxnorm_host[0];
  af::freeHost(maxnorm_host);
  return maxnorm;
}

// Minimum value 
double minval(const af::array& a){
  double *minval_host=NULL;
  minval_host = min(min(min(min(a,0),1),2),3).host<double>();
  double minval = minval_host[0];
  af::freeHost(minval_host);
  return minval;
}


//TODO check with c++14 (we used uncommented due to incompability with c++11 needed by cython)
////RK4 based on https://rosettacode.org/wiki/Runge-Kutta_method
//auto rk4(af::array f(double, af::array))
//{
//        return
//        [       f            ](double t,  af::array y, double dt ) -> af::array { return
//        [t,y,dt,f            ](                    af::array  dy1) -> af::array { return
//        [t,y,dt,f,dy1        ](                    af::array  dy2) -> af::array { return
//        [t,y,dt,f,dy1,dy2    ](                    af::array  dy3) -> af::array { return
//        [t,y,dt,f,dy1,dy2,dy3](                    af::array  dy4) -> af::array { return
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
//        auto eval_solution = [               ](double t          )->double{ return pow(t*t+4,2)/16                   ; } ;
//        auto find_error    = [eval_solution  ](double t, double y)->double{ return fabs(y-eval_solution(t))          ; } ;
//        auto is_whole      = [WHOLE_TOLERANCE](double t          )->bool  { return fabs(t-round(t)) < WHOLE_TOLERANCE; } ;
// 
//        auto dy = rk4( eval_diff_eqn ) ;
// 
//        double y = Y_START, t = T_START ;
// 
//        while(t <= TIME_MAXIMUM) {
//          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y,find_error(t,y)); }
//          y += dy(t,y,DT) ; t += DT;
//        }
//        return 0;
//}

