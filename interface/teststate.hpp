#ifndef TEST_STATE_H
#define TEST_STATE_H
#include "arrayfire.h"
#include "../src/mesh.hpp"
#include "../src/param.hpp"

class testState{
  public:
    testState (Mesh mesh_in, Param param_in, long int a);
    //testState (Mesh mesh_in, Param param_in, af::array a);
    af::array m;
    Mesh mesh;
    Param param;
    void printn0();
};
#endif

