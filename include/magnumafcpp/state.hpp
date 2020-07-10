#pragma once
#include "arrayfire.h"
#include "mesh.hpp"
#include <iostream>

namespace magnumafcpp {

class State {
  public:
    // Micromagnetic:
    State(Mesh mesh_in, double Ms, af::array m_in, bool verbose = true,
          bool mute_warning = false);
    State(Mesh mesh_in, af::array Ms_field, af::array m_in, bool verbose = true,
          bool mute_warning = false);
    // No Mesh:
    State(af::array m, double Ms, bool verbose = true,
          bool mute_warning = false);
    State(af::array m, af::array Ms_field, bool verbose = true,
          bool mute_warning = false);
    // Wrapping:
    State(Mesh mesh_in, double Ms, long int m_in, bool verbose = true,
          bool mute_warning = false);
    State(Mesh mesh_in, long int Ms_field_ptr, long int m_in,
          bool verbose = true, bool mute_warning = false);

    State operator+(const af::array&) const;
    void set_m(long int aptr); ///< For wrapping only: Setting member af::array
                               ///< m to values obtained from wrapped af.array
    long int get_m_addr();
    Mesh mesh{0, 0, 0, 0, 0, 0};
    double t{0.}; // time
    af::array m;  //!< magnetic field configuration
    double Ms{0}; //!< Saturation magnetization in [J/T/m^3]
    af::array
        Ms_field; //!< Non-homugenuous, mesh dependent saturation magnetization
                  //!< defined at every node in units of [J/T/m^3]. Is impicitly
                  //!< set and used when magnetization has values of norm 0.
    void set_Ms_field(long int afarray_ptr); // for wrapping only
    long int get_Ms_field();

    void set_Ms_field_if_m_minvalnorm_is_zero(const af::array& m,
                                              af::array& Ms_field);
    void check_m_norm(double tol = 1e-6);
    unsigned long long steps{0};
    void Normalize(); ///< normalize the magnetization to 1
    // af::array m_out;
    // long int get_m_addr(){m.lock(); return (long int) m.get();}

    void write_vti(std::string outputname);
    void _vti_writer_micro_boolean(std::string outputname);
    void _vti_writer_atom(std::string outputname);
    void _vti_reader(std::string inputname);

    double meani(const int i);
    void calc_mean_m(std::ostream& myfile);
    void calc_mean_m(std::ostream& myfile, double hzee);
    void calc_mean_m(std::ostream& myfile, const af::array& hzee);
    void calc_mean_m_steps(std::ostream& myfile, double hzee);
    void calc_mean_m_steps(std::ostream& myfile, const af::array& hzee);
    unsigned int get_n_cells_() { return n_cells_; };

    bool verbose{true};
    bool mute_warning{false};
    bool afsync{false};

  private:
    ///< Number of cells with Ms != 0
    unsigned int n_cells_{0};
    ///< Boolean array of type b8 and size [x, y, z, 1] indicating whether the
    ///< respective cell is considered in mean value calulation (==1) or not
    ///< (==0)
    af::array evaluate_mean_;
    ///< Number of cells with for which evaluate_mean_ is 1
    unsigned int evaluate_mean_is_1_{0};

    // void check_discretization();
    // void check_nonequispaced_discretization();
};
} // namespace magnumafcpp
