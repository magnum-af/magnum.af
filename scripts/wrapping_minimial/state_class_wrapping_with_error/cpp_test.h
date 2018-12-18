#ifndef TEST_H
#define TEST_H

#include<arrayfire.h>

class Test {
public:
    af::array m;
    Test();
    ~Test();
    void init_m(long int aptr);
    void print_m();
    long int get_m();
};
#endif
