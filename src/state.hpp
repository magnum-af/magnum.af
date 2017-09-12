#ifndef STATE_H
#define STATE_H
#include "arrayfire.h"
#include "mesh.hpp"
#include "param.hpp"

class State{
  public:
    State (Mesh mesh_in, Param param_in, af::array m_in);
    Mesh mesh;
    Param param;
    double t{0.};//time
    af::array m;
    int steps{0};
    //State& operator+(af::array& rhs);
};

#endif
