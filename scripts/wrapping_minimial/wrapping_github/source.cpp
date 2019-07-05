#include "source.hpp"

using namespace magnumaf;


Wrap:: Wrap() {}
Wrap::~Wrap() {}

void Wrap::c_py_to_cpp(long int addr){
    A = *( new af::array( *(void**) addr ));
}

void Wrap::c_calc_B(){
    af::print("cpp: A= ", A);
    B=2*A;
    af::print("cpp: B= ", B);
}

long int Wrap::c_cpp_to_py(){
    B.lock();
    return (long int) B.get();
}
