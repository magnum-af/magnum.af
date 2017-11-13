//#pragma GCC clib library_with_useful_functions
#ifndef STATE_H
#define STATE_H
#include "arrayfire.h"
#include "mesh.hpp"
#include "param.hpp"

class State{
  public:
    State (Mesh mesh_in, Param param_in, af::array m_in);
    State (Mesh mesh_in, Param param_in, long int aptr);
    Mesh mesh;
    Param param;
    double t{0.};//time
    af::array m;
    int steps{0};
    void print_m();//TODO Replace with return array


    //State& operator+(af::array& rhs);
};

#endif
