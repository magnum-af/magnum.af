#include "source.hpp"

Wrap:: Wrap() {}
Wrap::~Wrap() {}

void Wrap::c_py_to_cpp(long int addr){
    A = *( new af::array( *(void**) addr ));
}

void Wrap::c_calc_B(){
    af::print("cpp: A= ", A);
    B=2*A;
    af::print("cpp: B= ", B);
    //with memory problems try
    //B.lock();
}

long int Wrap::c_cpp_to_py(){
    return (long int) B.get();
}

