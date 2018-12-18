#include "cpp_test.h"
#include <stdio.h>

Test::Test() { }

Test::~Test() 
{
    printf("Cleanup!\n");
}

void Test::init_m(long int aptr){
    void **a = (void **)aptr;
    af::array *A = new af::array(*a);
    m = *A;
}

void Test::print_m(){
    af::print("m=", m);
}

long int Test::get_m() {
    return (long int) this->m.get();
}
