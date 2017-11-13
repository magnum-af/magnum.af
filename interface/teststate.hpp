//#include "arrayfire.h"
#include "../src/mesh.hpp"
#include "../src/param.hpp"

class testState{
  public:
    testState (Mesh mesh_in, Param param_in);
    Mesh mesh;
    Param param;
};

