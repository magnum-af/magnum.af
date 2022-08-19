#pragma once
#include "arrayfire.h"
#include "mesh.hpp"
#include "util/util.hpp"
#include <array>
#include <iostream>
#include <limits>

namespace magnumaf {

class State {
  public:
    // Micromagnetic:
    State(Mesh mesh_in, double Ms, const af::array& m_in, bool verbose = true, bool mute_warning = false);
    State(Mesh mesh_in, af::array Ms_field, const af::array& m_in, bool verbose = true, bool mute_warning = false);
    // No Mesh:
    State(const af::array& m, double Ms, bool verbose = true, bool mute_warning = false);
    State(const af::array& m, af::array Ms_field, bool verbose = true, bool mute_warning = false);
    // Wrapping:
    State(Mesh mesh_in, double Ms, long int m_in, bool verbose = true, bool mute_warning = false);
    State(Mesh mesh_in, long int Ms_field_ptr, long int m_in, bool verbose = true, bool mute_warning = false);

    State operator+(const af::array&) const;
    /// Writes <mx> <my> <mz> to stream, separated by tabs.
    friend std::ostream& operator<<(std::ostream& os, const State&);

    void set_m(long int aptr); ///< For wrapping only: Setting member af::array
                               ///< m to values obtained from wrapped af.array
    long int get_m_addr() const;
    Mesh mesh{0, 0, 0, 0, 0, 0};
    double t{0.};                                        // time
    af::array m;                                         //!< magnetic field configuration
    double Ms{std::numeric_limits<double>::quiet_NaN()}; //!< Saturation magnetization in [J/T/m^3] or [A/m]
    af::array Ms_field;                                  //!< Inhomogeneous, mesh dependent saturation magnetization
                                                         //!< defined at every node in units of [J/T/m^3]. Is impicitly
                                                         //!< set and used when magnetization has values of norm 0.

    af::array get_Ms_field_in_vec_dims() const { return af::tile(Ms_field, 1, 1, 1, 3); }
    af::array get_Ms_as_field()
        const; //!< return Ms as af::array. If Ms_field is not set, scalar Ms is tiled to dims [nx, ny, nz, 1].
    af::array get_Ms_as_field_in_vector_dims() const; //!< return Ms tiled to dims [nx, ny, nz, 3].
    void set_Ms_field(long int afarray_ptr);          // for wrapping only
    long int wrapping_get_Ms_field() const;

    void set_Ms_field_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms_field);
    void check_m_norm(double tol = 1e-6) const;
    unsigned long long steps{0};
    void Normalize(); ///< normalize the magnetization to 1
    // af::array m_out;
    // long int get_m_addr(){m.lock(); return (long int) m.get();}

    void write_vti(const std::string& outputname) const;
    void _vti_writer_atom(std::string outputname) const;
    void _vti_reader(const std::string& inputname);

    /// Get the i'th component of <m>: 0 == mx, 1 == my, 2 == mz
    double meani(const int i) const;
    /// Returns {<mx>, <my>, <mz>}, the average magnetization in each spacial direction
    std::array<double, 3> mean_m() const;
    af::array mean_m_as_afarray() const;
    /// Returns <mx>, i.e. average magnetisation in x-direction
    double mean_mx() const { return mean_m()[0]; }
    /// Returns <my>, i.e. average magnetisation in y-direction
    double mean_my() const { return mean_m()[1]; }
    /// Returns <mz>, i.e. average magnetisation in z-direction
    double mean_mz() const { return mean_m()[2]; }

    /// Spacial Average Magnetic Moment per volume <M> in [J/T/m^3] or [A/m],
    /// defined as ( /sum_i m_i Ms_i ) / N, where N is no of cells where |m| != 0
    auto mean_M_as_afarray() const -> af::array;
    long int wrapping_mean_M_as_afarray() const { return util::pywrap::send_copy_to_py(mean_M_as_afarray()); }
    auto mean_M() const -> std::array<double, 3>;

    bool verbose{true};
    bool mute_warning{false};

  private:
};
} // namespace magnumaf
