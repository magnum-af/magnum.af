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
    void set_Ms_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms);
    void check_discretization();
    int steps{0};
    long int get_m_addr(){return (long int) m.get();}
    //af::array m_out;
    //long int get_m_addr(){m.lock(); return (long int) m.get();}

    void _vti_writer_micro(std::string outputname);
    void _vti_writer_atom (std::string outputname);
    void _vti_reader(std::string inputname);

    void _vtr_writer(std::string outputname);
    void _vtr_reader(std::string inputname);
    double meani(const int i);
    void calc_mean_m( std::ostream& myfile);
    void calc_mean_m( std::ostream& myfile, const long int n_cells ); // n_cells is number of cells with non_zero_Ms
    void calc_mean_m( std::ostream& myfile, const long int n_cells, double hzee);
    void calc_mean_m( std::ostream& myfile, const long int n_cells, const af::array& hzee);
};

#endif
