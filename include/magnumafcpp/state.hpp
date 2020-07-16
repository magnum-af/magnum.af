#pragma once
#include "arrayfire.h"
#include "mesh.hpp"
#include <array>
#include <iostream>

namespace magnumafcpp {

class State {
  public:
    // Micromagnetic:
    State(Mesh mesh_in, double Ms, af::array m_in, bool verbose = true, bool mute_warning = false);
    State(Mesh mesh_in, af::array Ms_field, af::array m_in, bool verbose = true, bool mute_warning = false);
    // No Mesh:
    State(af::array m, double Ms, bool verbose = true, bool mute_warning = false);
    State(af::array m, af::array Ms_field, bool verbose = true, bool mute_warning = false);
    // Wrapping:
    State(Mesh mesh_in, double Ms, long int m_in, bool verbose = true, bool mute_warning = false);
    State(Mesh mesh_in, long int Ms_field_ptr, long int m_in, bool verbose = true, bool mute_warning = false);

    State operator+(const af::array&) const;
    /// Writes <mx> <my> <mz> to stream, separated by tabs.
    friend std::ostream& operator<<(std::ostream& os, const State&);

    void set_m(long int aptr); ///< For wrapping only: Setting member af::array
                               ///< m to values obtained from wrapped af.array
    long int get_m_addr();
    Mesh mesh{0, 0, 0, 0, 0, 0};
    double t{0.};                            // time
    af::array m;                             //!< magnetic field configuration
    double Ms{0};                            //!< Saturation magnetization in [J/T/m^3]
    af::array Ms_field;                      //!< Non-homugenuous, mesh dependent saturation magnetization
                                             //!< defined at every node in units of [J/T/m^3]. Is impicitly
                                             //!< set and used when magnetization has values of norm 0.
    void set_Ms_field(long int afarray_ptr); // for wrapping only
    long int get_Ms_field();

    void set_Ms_field_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms_field);
    void check_m_norm(double tol = 1e-6);
    unsigned long long steps{0};
    void Normalize(); ///< normalize the magnetization to 1
    // af::array m_out;
    // long int get_m_addr(){m.lock(); return (long int) m.get();}

    void write_vti(std::string outputname);
    void _vti_writer_atom(std::string outputname);
    void _vti_reader(std::string inputname);

    /// Get the i'th component of <m>: 0 == mx, 1 == my, 2 == mz
    double meani(const int i);
    /// Returns {<mx>, <my>, <mz>}, the average magnetization in each spacial direction
    std::array<double, 3> mean_m() const;
    af::array mean_m_as_afarray() const;
    /// Returns <mx>, i.e. average magnetisation in x-direction
    double mean_mx() { return mean_m()[0]; }
    /// Returns <my>, i.e. average magnetisation in y-direction
    double mean_my() { return mean_m()[1]; }
    /// Returns <mz>, i.e. average magnetisation in z-direction
    double mean_mz() { return mean_m()[2]; }

    unsigned int get_n_cells_() { return n_cells_; };

    bool verbose{true};
    bool mute_warning{false};
    bool afsync{false};

  private:
    ///< Number of cells with Ms != 0
    unsigned int n_cells_{0};
};
} // namespace magnumafcpp
