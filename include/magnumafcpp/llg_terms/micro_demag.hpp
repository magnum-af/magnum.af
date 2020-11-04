#pragma once
#include "../state.hpp"
#include "arrayfire.h"
#include "integrator_term_mesh_base.hpp"

namespace magnumafcpp {

namespace newell {
double Nxx(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz);
double Nxy(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz);
} // namespace newell

class DemagField : public IntegratorTermMeshBase {
  public:
    DemagField(Mesh, bool verbose = false, bool caching = true, unsigned nthreads = 0);

    af::array h(const State& state);

    double get_cpu_time() { return cpu_time; }

    ///< Get copy of array storing the Fourier transfrom of the demag tensor.
    af::array get_Nfft() const { return Nfft; }
    // For wrapping
    void print_Nfft() const;

  private:
    af::array Nfft;
    double cpu_time{0.};
};
} // namespace magnumafcpp
