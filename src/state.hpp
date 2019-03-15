#ifndef STATE_H
#define STATE_H
#include "arrayfire.h"
#include "mesh.hpp"
#include "material.hpp"
#include "vtk_IO.hpp"
#include "func.hpp"
#include "misc.hpp"


class State{
  public:
    State (Mesh mesh_in, Material param_in, af::array m_in);
    State (Mesh mesh_in, Material param_in, af::array m_in, af::array evaluate_mean);
    State (Mesh mesh_in, Material param_in, long int aptr);
    State (Mesh mesh_in, Material param_in, long int aptr, long int evaluate_mean_ptr);
    ~State(){};
    void set_m(long int aptr); ///< For wrapping only: Setting member af::array m to values obtained from wrapped af.array
    long int get_m_addr();
    Mesh mesh;
    Material material;
    double t{0.};//time
    af::array m;
    af::array Ms; // Saturation magnetization. Is impicitly set and used when magnetization has values of norm 0.
    void set_micro_Ms_field(long int aptr);
    long int get_micro_Ms_field();

    af::array micro_A_field; //!< Spacially varying exchange energy in [J/m] defined at each node.//TODO move to material or exchange
    void set_micro_A_field(long int aptr); ///< For wrapping only: Setting af::array micro_A_field.
    long int get_micro_A_field();

    af::array micro_Ku1_field; //!< Spacially varying anisotropy energy in [J/m^3] defined at each node.//TODO move to material or exchange
    void set_micro_Ku1_field(long int aptr); ///< For wrapping only: Setting af::array micro_Ku1_field.
    long int get_micro_Ku1_field();

    void set_Ms_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms);
    void check_discretization();
    void check_m_norm(double tol = 1e-6);
    unsigned long long steps{0};
    void Normalize(); ///< normalize the magnetization to 1
    //af::array m_out;
    //long int get_m_addr(){m.lock(); return (long int) m.get();}

    void _vti_writer_micro(std::string outputname);
    void _vti_writer_micro_boolean(std::string outputname);
    void _vti_writer_atom (std::string outputname);
    void _vti_reader(std::string inputname);

    void _vtr_writer(std::string outputname);
    void _vtr_reader(std::string inputname);
    double meani(const int i);
    void calc_mean_m( std::ostream& myfile);
    void calc_mean_m( std::ostream& myfile, double hzee);
    void calc_mean_m( std::ostream& myfile, const af::array& hzee);
    void calc_mean_m_steps( std::ostream& myfile, double hzee);
    void calc_mean_m_steps( std::ostream& myfile, const af::array& hzee);
    unsigned int get_n_cells_(){return n_cells_;};

    bool verbose_{true};
  private:
    ///< Number of cells with Ms != 0
    unsigned int n_cells_{0};
    ///< Boolean array of type b8 and size [x,y,z,1] indicating whether the respective cell is considered in mean value calulation (==1) or not (==0)
    af::array evaluate_mean_;
    ///< Number of cells with for which evaluate_mean_ is 1
    unsigned int evaluate_mean_is_1_{0};

};

#endif
