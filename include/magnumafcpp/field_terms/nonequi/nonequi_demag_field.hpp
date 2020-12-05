#pragma once
#include "state.hpp"
#include "arrayfire.h"
#include "nonequi_term_base.hpp"

namespace magnumafcpp {

class NonequiDemagField : public NonequiTermBase {
  public:
    NonequiDemagField(NonequiMesh nonequimesh, bool verbose = true, bool caching = false, unsigned nthreads = 0);
    const af::array Nfft; //!< Array storing the Fourier transfrom of the demag tensor.

    af::array h(const State& state) const override;           // Field contribution
    void print_Nfft();                         // For wrapping

  private:
    double cpu_time{0.};
};

// expose for testing:
namespace newell_nonequi {
double Nxx(const double x, const double y, const double z, const double dx, const double dy, const double dz,
           const double dX, const double dY, const double dZ);
double Nxy(const double x, const double y, const double z, const double dx, const double dy, const double dz,
           const double dX, const double dY, const double dZ);
double nonequi_index_distance(const std::vector<double> spacings, const unsigned i, const unsigned j,
                              const bool verbose = true);

} // namespace newell_nonequi
} // namespace magnumafcpp
