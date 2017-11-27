#ifndef STATE_H
#define STATE_H
#include "arrayfire.h"
#include "mesh.hpp"
#include "param.hpp"

class State{
  public:
    State (Mesh mesh_in, Param param_in, af::array m_in);
    State (Mesh mesh_in, Param param_in, long int aptr);
    ~State(){};
    Mesh mesh;
    Param param;
    double t{0.};//time
    af::array m;
    int steps{0};
    //long int get_m_addr(){return (long int) m.get();}
    long int get_m_addr(){m.lock(); return (long int) m.get();}


    //State& operator+(af::array& rhs);
};

#endif
