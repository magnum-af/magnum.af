#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumafcpp {

namespace newell {
double Nxx(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz);
double Nxy(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz);
} // namespace newell

class DemagField : public MicroTerm {
  public:
    DemagField(Mesh, bool verbose = false, bool caching = true, unsigned nthreads = 0);

    virtual af::array h(const State& state) const override;

    ///< Get copy of array storing the Fourier transfrom of the demag tensor.
    af::array get_Nfft() const { return Nfft; }
    // For wrapping
    void print_Nfft() const;

  private:
    mutable af::array Nfft; // mutable for c64-c32 conversion
    double cpu_time{0.};
};
} // namespace magnumafcpp
