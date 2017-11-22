#ifndef SOURCE_H
#define SOURCE_H
#include<arrayfire.h>

class Wrap {
public:
    Wrap();
    ~Wrap();
    
    void c_py_to_cpp(long int a_addr);
    long int c_cpp_to_py();
    void c_calc_B();

    af::array A;
    af::array B;
};
#endif
