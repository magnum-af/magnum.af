#ifndef STATE_H
#define STATE_H
#include "arrayfire.h"
#include "mesh.hpp"
#include "param.hpp"
#include "vtk_IO.hpp"


class State{
  public:
    State (Mesh mesh_in, Param param_in, af::array m_in);
    State (Mesh mesh_in, Param param_in, long int aptr);
    ~State(){};
    Mesh mesh;
    Param param;
    double t{0.};//time
    af::array m;
    af::array Ms; // Saturation magnetization
    af::array Inverse_Ms; // Inverse Saturation magnetization Inverse_Ms=1/Ms: to avoid division by Ms with zero
    int steps{0};
    long int get_m_addr(){return (long int) m.get();}
    //long int get_m_addr(){m.lock(); return (long int) m.get();}

    void _vti_writer_micro(std::string outputname);
    void _vti_writer_atom (std::string outputname);
    void _vti_reader(std::string inputname);

    void _vtr_writer(std::string outputname);
    void _vtr_reader(std::string inputname);
};

#endif
