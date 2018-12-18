#ifndef TEST_H
#define TEST_H

#include<arrayfire.h>

class Test {
  public:
    Test();
    ~Test();
    void initialize_m(long int aptr);
    void manipulate_m();
    long int get_m();
  private:
    af::array m;
};
#endif
