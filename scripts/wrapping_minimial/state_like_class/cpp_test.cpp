#include <stdio.h>
#include "cpp_test.h"

using namespace magnumaf;


Test::Test() { }

Test::~Test()
{
    printf("Cleanup!\n");
}

void Test::initialize_m(long int aptr){
    void **a = (void **)aptr;
    af::array *A = new af::array(*a);
    this->m = *A;
}

void Test::manipulate_m(){
    this->m *= 2.;
}

long int Test::get_m() {
    af::array* a = new af::array(this->m);
    return (long int) a->get();
}
