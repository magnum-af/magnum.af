#pragma once
#include "arrayfire.h"
#include "field_terms/micro/micro_term.hpp"
#include "state.hpp"

namespace magnumaf {

// Exporting functions for unit tests
namespace newell {
double Nxx(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz);
double Nxy(const int ix, const int iy, const int iz, const double dx, const double dy, const double dz);
} // namespace newell

class DemagField : public MicroTerm {
  public:
    DemagField(Mesh, bool verbose = false, bool caching = true, unsigned nthreads = 0);

    ///< Get copy of array storing the Fourier transfrom of the demag tensor.
    af::array get_Nfft() const { return this->Nfft; }
    // For wrapping
    long int get_Nfft_ptr() const { return util::pywrap::send_copy_to_py(this->Nfft); }
    void print_Nfft() const { af::print("Nfft=", Nfft); }

  private:
    mutable af::array Nfft; // mutable for c64-c32 conversion
    virtual af::array impl_H_in_Apm(const State& state) const override;
};
} // namespace magnumaf
