#include "gradient_decent.hpp"
#include<math.h>
#include<iostream>

namespace magnumafcpp{
//magnumafcpp::GradientDecent::GradientDecent(std::function<double (double)> f, Precision precision){
//    std::cout << precision.get() << std::endl;
//}

double GradientDecent::calc_df(double x_current){
    double f_x_current = f(x_current);
    if(verbose.get()) std::cout << "calc_df: x_current=" << x_current << ", f_x_current=" << f_x_current << std::endl;
    return ( f(x_current + epsilon.get()) - f_x_current )/epsilon.get();
}

std::pair<double, double> GradientDecent::minimize(){
    //double df = ( f(x_start_val.get()) - f(x_start_val.get() + epsilon.get()) )/epsilon.get();
    double x_min = 1e-10;
    double x_next = x_start_val.get();
    for (int i = 0; i < maxiters.get(); i ++){
        double x_current = x_next;
        double df = calc_df(x_current);
        x_next = x_current - gamma.get() * df;
        //x_next = x_current - gamma.get() * calc_df(x_current);

        if (x_next <= x_min){
            printf("Warning, x_nex <= x_min, setting x_next to x_min\n");
            x_next = x_min;
        }
        double step = x_next - x_current;
        if(verbose.get()) std::cout << "minimi: x_current=" << x_current << ", x_next= " << x_next << ", step=" << step << ", df=" << df << std::endl;
        if (std::fabs(step) < precision.get()) break;
    }

    return std::pair<double, double> (x_next, f(x_next));
}

}// namespace magnumafcpp
